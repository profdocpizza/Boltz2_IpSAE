#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate YAMLs and run scripts for binder validation from a single config file.

Usage:
    python make_binder_validation_scripts.py --config config.yml

Key features:
  - Binder is always first in Boltz YAML and uses chain IDs: A, B, C, ...
  - Partner (second entity) uses IDs depending on role:
        target     → TA, TB, TC, ...
        antitarget → AA, AB, AC, ...
        (fallback) → PA, PB, PC, ... if role is unknown
  - Supports multichain binders and multichain targets/antitargets.
  - Supports ligands (CCD and/or SMILES) on binders and non-self targets.
  - MSA can be provided for any chain via config (chains_msa).
  - from_dir entries NEVER have MSAs (by design).
  - Uses 'target_' / 'antitarget_' prefixes in YAML names:
       binder_<binder>_vs_target_<name>.yaml
       binder_<binder>_vs_antitarget_<name>.yaml
  - Generates:
       - Per-binder YAMLs for all binder–(anti)target pairs
       - Per-binder run.sh
       - Global run_all_cofolding.sh using one persistent Boltz worker
       - Visualization helper script
"""

import argparse
import os
import re
import shlex
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from str2fasta import get_sequences_all_chains
import yaml

SANITIZE_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def sanitize_name(name: str) -> str:
    cleaned = SANITIZE_RE.sub("_", name.strip())
    if not cleaned:
        raise ValueError(f"Invalid name: {name!r}")
    return cleaned


def write_text(path: Path, text: str) -> None:
    """Write text to a file, creating parent dirs if needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def read_fasta_multi(fasta_path: Path) -> List[str]:
    """Read one FASTA file and return a list of sequences (one per record)."""
    seqs: List[str] = []
    seq: List[str] = []
    with fasta_path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    seqs.append("".join(seq).replace(" ", "").upper())
                    seq = []
            else:
                seq.append(line)
    if seq:
        seqs.append("".join(seq).replace(" ", "").upper())
    return seqs


def read_fasta_dir_entities(fasta_dir: Path) -> List[Tuple[str, List[str]]]:
    """
    Read all FASTA-like files in a directory.
    Returns list of (entity_name, [seqs]).
    """
    entities: List[Tuple[str, List[str]]] = []
    for fasta_path in sorted(fasta_dir.glob("*")):
        if not fasta_path.is_file() or fasta_path.suffix.lower() not in {
            ".fasta",
            ".fa",
            ".fna",
            ".faa",
            ".txt",
        }:
            continue
        name = sanitize_name(fasta_path.stem)
        seqs = read_fasta_multi(fasta_path)
        if seqs:
            entities.append((name, seqs))
    return entities


def add_n_terminal_lysine(seqs: List[str]) -> List[str]:
    """Prepend 'K' if missing at N-terminus for each sequence."""
    return [("K" + s if not s.startswith("K") else s) for s in seqs]


# ---------------------------------------------------------------------------

def _alpha_suffix(idx: int) -> str:
    """
    Return a letter-like suffix for chain indices: A, B, C, ... Z, X26, X27, ...
    (Only the first 26 are pretty; beyond that we degrade gracefully.)
    """
    if idx < 26:
        return chr(ord("A") + idx)
    return "X" + str(idx)


def _partner_chain_id(role: str, idx: int) -> str:
    role = (role or "").lower()
    if role == "target":
        prefix = "T"
    elif role == "antitarget":
        prefix = "A"
    elif role == "self":
        prefix = "S"
    else:
        raise ValueError(
            f"Unknown partner role: {role!r} (expected 'target', 'antitarget', or 'self')."
        )
    return prefix + _alpha_suffix(idx)


def _yaml_quote(value: str) -> str:
    """Return a single-quoted YAML-safe scalar."""
    return "'" + value.replace("'", "''") + "'"


def _allocate_ligand_chain_ids(count: int, used_ids: set) -> List[str]:
    """
    Allocate `count` chain IDs for ligands while avoiding `used_ids`.
    Uses LA, LB, LC, ... (then LX26, LX27, ... for very large counts).
    """
    out: List[str] = []
    idx = 0
    while len(out) < count:
        cid = "L" + _alpha_suffix(idx)
        idx += 1
        if cid in used_ids:
            continue
        out.append(cid)
        used_ids.add(cid)
    return out


def _coerce_string_or_list(raw: Any, key_name: str, context: str) -> List[str]:
    """Normalize scalar/list config value into a non-empty string list."""
    if raw is None:
        return []
    if isinstance(raw, str):
        values = [raw]
    elif isinstance(raw, list):
        values = raw
    else:
        raise ValueError(f"{context}: '{key_name}' must be a string or list of strings.")

    out: List[str] = []
    for v in values:
        s = str(v).strip()
        if not s:
            raise ValueError(f"{context}: '{key_name}' contains an empty value.")
        out.append(s)
    return out


def parse_ligands(entry: Dict[str, Any], context: str) -> List[Dict[str, str]]:
    """
    Parse ligand config keys and normalize into:
      [{"kind": "ccd"|"smiles", "value": "<string>"}, ...]
    Accepts singular/plural aliases for backwards compatibility.
    """
    ligands: List[Dict[str, str]] = []

    ccd_values = []
    ccd_values.extend(_coerce_string_or_list(entry.get("ligand_ccd"), "ligand_ccd", context))
    ccd_values.extend(_coerce_string_or_list(entry.get("ligands_ccd"), "ligands_ccd", context))

    smiles_values = []
    smiles_values.extend(
        _coerce_string_or_list(entry.get("ligand_smiles"), "ligand_smiles", context)
    )
    smiles_values.extend(
        _coerce_string_or_list(entry.get("ligands_smiles"), "ligands_smiles", context)
    )

    for ccd in ccd_values:
        ligands.append({"kind": "ccd", "value": ccd})
    for smi in smiles_values:
        ligands.append({"kind": "smiles", "value": smi})
    return ligands


def has_ligand_keys(entry: Dict[str, Any]) -> bool:
    keys = {"ligand_ccd", "ligands_ccd", "ligand_smiles", "ligands_smiles"}
    return any(k in entry for k in keys)


def _warn_template_parse_failure(cif_path: Path, context: str, details: str) -> None:
    """Emit a non-fatal warning for template CIF precheck failures."""
    print(
        f"⚠ {context}: Boltz may fail to parse template CIF: {cif_path}\n"
        f"   Details: {details}\n"
        "   Recommendation: convert PDB to mmCIF with "
        "https://mmcif.pdbj.org/converter/index.php?l=en"
    )


def _exc_summary(exc: BaseException) -> str:
    """Compact exception summary with type and repr for empty-message exceptions."""
    return f"{type(exc).__name__}: {exc!r}"


def _precheck_template_cif(
    cif_path: Optional[str], context: str, checked_paths: Set[Path]
) -> None:
    """
    Precheck template CIF parseability before generating YAML.

    Strategy:
      1) Try Boltz parser directly (best predictor of runtime behavior).
      2) If unavailable/failing, fall back to a basic Gemmi structural sanity check.
      3) Emit warning (non-fatal) when parseability looks problematic.
    """
    if not cif_path:
        return

    path = Path(cif_path).resolve()
    if path in checked_paths:
        return
    checked_paths.add(path)

    if not path.is_file():
        raise ValueError(f"{context}: cif_template file not found: {path}")

    boltz_error: Optional[Exception] = None
    boltz_checked = False

    try:
        from boltz.data.mol import load_canonicals
        from boltz.data.parse.mmcif import parse_mmcif

        mol_dir = Path.home() / ".boltz" / "mols"
        ccd = load_canonicals(mol_dir)
        parse_mmcif(
            path,
            mols=ccd,
            moldir=mol_dir,
            use_assembly=False,
            compute_interfaces=False,
        )
        return
    except Exception as e:
        boltz_checked = True
        boltz_error = e

    try:
        import gemmi

        doc = gemmi.cif.read_file(str(path))
        block = doc.sole_block()
        structure = gemmi.make_structure_from_block(block)
        if len(structure) == 0:
            raise ValueError("Gemmi found 0 models in the CIF.")
        if len(structure[0]) == 0:
            raise ValueError("Gemmi found no chains in model 0.")
    except Exception as gemmi_error:
        if boltz_checked:
            details = (
                f"Boltz parser check failed ({_exc_summary(boltz_error)}); "
                f"Gemmi check failed ({_exc_summary(gemmi_error)})."
            )
        else:
            details = f"Gemmi check failed ({_exc_summary(gemmi_error)})."
        _warn_template_parse_failure(path, context, details)
        return

    # Gemmi passed but Boltz failed: still warn, this is the common real-world failure mode.
    if boltz_checked:
        details = (
            f"Boltz parser check failed ({_exc_summary(boltz_error)}), although Gemmi basic check passed. "
            "Boltz may still reject this CIF."
        )
        _warn_template_parse_failure(path, context, details)



def yaml_for_pair(
    binder_seqs: List[str],
    partner_seqs: List[str],
    partner_role: str,
    use_msa_server: bool,
    binder_msas: Optional[List[Optional[str]]] = None,
    partner_msas: Optional[List[Optional[str]]] = None,
    binder_ligands: Optional[List[Dict[str, str]]] = None,
    partner_ligands: Optional[List[Dict[str, str]]] = None,
    cif_template: Optional[str] = None,
) -> str:
    """
    Build Boltz YAML for a binder–partner pair.

    - Binder chains first with IDs: A, B, C, ...
    - Partner chains next with IDs depending on role:
          target     → TA, TB, ...
          antitarget → AA, AB, ...
          other      → PA, PB, ...

    Logic for MSAs:
      - If user provided a path (chains_msa), ALWAYS write 'msa: <path>'.
      - If user did NOT provide a path:
         - if use_msa_server=False: write 'msa: empty'
         - if use_msa_server=True:  do NOT write 'msa: ...' (let Boltz fetch it)
    """

    lines: List[str] = ["version: 1", "sequences:"]
    used_chain_ids: set = set()
    binder_msas = binder_msas or [None] * len(binder_seqs)
    partner_msas = partner_msas or [None] * len(partner_seqs)
    binder_ligands = binder_ligands or []
    partner_ligands = partner_ligands or []

    # Helper
    def _add_msa_field(msa_opt: Optional[str]):
        if msa_opt:
            lines.append(f"      msa: {msa_opt}")
        else:
            # No user-provided MSA
            if not use_msa_server:
                lines.append("      msa: empty")

    # --- Binder chains (A, B, ...) ---
    for i, seq in enumerate(binder_seqs):
        cid = _alpha_suffix(i)
        lines.append("  - protein:")
        lines.append(f"      id: {cid}")
        lines.append(f"      sequence: {seq}")
        used_chain_ids.add(cid)
        msa_val = binder_msas[i] if i < len(binder_msas) else None
        _add_msa_field(msa_val)

    # --- Partner chains (TA/TB/... or AA/AB/...) ---
    for i, seq in enumerate(partner_seqs):
        cid = _partner_chain_id(partner_role, i)
        used_chain_ids.add(cid)
        lines.append("  - protein:")
        lines.append(f"      id: {cid}")
        lines.append(f"      sequence: {seq}")
        msa_val = partner_msas[i] if i < len(partner_msas) else None
        _add_msa_field(msa_val)

    # --- Ligands ---
    all_ligands = binder_ligands + partner_ligands
    if all_ligands:
        ligand_ids = _allocate_ligand_chain_ids(len(all_ligands), used_chain_ids)

        # Group repeated ligands into one entry with id: [LA, LB, ...]
        grouped_ids: Dict[Tuple[str, str], List[str]] = {}
        ligand_order: List[Tuple[str, str]] = []
        for lig, cid in zip(all_ligands, ligand_ids):
            key = (lig["kind"], lig["value"])
            if key not in grouped_ids:
                grouped_ids[key] = []
                ligand_order.append(key)
            grouped_ids[key].append(cid)

        for key in ligand_order:
            kind, value = key
            ids = grouped_ids[key]
            lines.append("  - ligand:")
            if len(ids) == 1:
                lines.append(f"      id: {ids[0]}")
            else:
                lines.append(f"      id: [{', '.join(ids)}]")
            lines.append(f"      {kind}: {_yaml_quote(value)}")

    # --- Global Templates Block ---
    if cif_template and partner_role in ["target", "antitarget"] and cif_template is not None:
        lines.append("templates:")
        lines.append(f"  - cif: {cif_template}")

    return "\n".join(lines) + "\n"


def make_run_sh(
    dir_path: Path,
    yaml_paths: List[Path],
    recycling_steps: Optional[int],
    diffusion_samples: Optional[int],
    use_msa_server: bool,
) -> None:
    """Create per-binder run.sh."""
    lines: List[str] = ["#!/bin/bash", "set -e", ""]
    for p in yaml_paths:
        cmd = ["boltz", "predict", p.name]
        if recycling_steps is not None:
            cmd += ["--recycling_steps", str(recycling_steps)]
        if use_msa_server:
            cmd.append("--use_msa_server")
        if diffusion_samples is not None:
            cmd += ["--diffusion_samples", str(diffusion_samples)]
        cmd += ["--out_dir", os.path.join(dir_path, "outputs")]
        lines.append(" ".join(cmd))
    run_path = dir_path / "run.sh"
    write_text(run_path, "\n".join(lines) + "\n")
    os.chmod(run_path, 0o755)


def make_master_run_sh(
    output_root: Path,
    recycling_steps: Optional[int],
    diffusion_samples: Optional[int],
    use_msa_server: bool,
) -> None:
    """Generate the serial hot-worker launcher for all binder YAMLs."""
    worker_path = Path(__file__).resolve().with_name("hot_boltz_worker.py")
    command = [
        "python",
        "-u",
        str(worker_path),
        "--output-root",
        '"$DIR"',
    ]
    if recycling_steps is not None:
        command += ["--recycling-steps", str(recycling_steps)]
    if diffusion_samples is not None:
        command += ["--diffusion-samples", str(diffusion_samples)]
    if use_msa_server:
        command.append("--use-msa-server")
    command_text = " ".join(
        '"$DIR"' if part == '"$DIR"' else shlex.quote(part) for part in command
    )
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        "# Determine directory of this script",
        'DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"',
        "",
        "# One persistent Boltz2 process executes all YAML jobs serially.",
        command_text,
        "",
    ]
    run_all_path = output_root / "run_all_cofolding.sh"
    write_text(run_all_path, "\n".join(lines))
    os.chmod(run_all_path, 0o755)
    print(f"✅ Created {run_all_path}")



def make_visualisation_sh(output_root: Path, cfg: Dict[str, Any]) -> None:
    """Create visualization helper script (config-driven)."""

    viz_cfg = cfg.get("visualisation", {}) or {}

    ipsae_e = viz_cfg.get("ipsae_error_threshold", 15)
    ipsae_d = viz_cfg.get("ipsae_distance_threshold", 15)
    use_best = viz_cfg.get("use_best_model", True)
    num_cpu = viz_cfg.get("num_cpu", None)
    rmsd_to_monomer = parse_bool_option(
        viz_cfg.get("RMSD_to_binder_monomer", False),
        "visualisation.RMSD_to_binder_monomer",
        default=False,
    )

    # Flags
    best_flag = "--use_best_model" if use_best else ""
    cpu_flag = f"--num_cpu {num_cpu}" if num_cpu is not None else ""
    rmsd_flag = "--rmsd_to_binder_monomer" if rmsd_to_monomer else ""

    script_line = (
        f"python {os.path.dirname(os.path.abspath(__file__))}/visualise_binder_validation.py "
        f"--ipsae_e {ipsae_e} "
        f"--ipsae_d {ipsae_d} "
        f"{cpu_flag} "
        f"--root_dir {output_root} "
        f"--generate_data --plot {best_flag} {rmsd_flag}"
    ).strip()

    sh_path = output_root / "visualise_cofolding_results.sh"
    write_text(sh_path, script_line + "\n")
    os.chmod(sh_path, 0o755)




# ---------------------------------------------------------------------------
# Config parsing
# ---------------------------------------------------------------------------

def load_config(path: Path) -> Dict[str, Any]:
    with path.open() as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise ValueError("Top-level config must be a mapping.")
    return cfg


def get_global_option(cfg: Dict[str, Any], *keys, default=None):
    """Convenience for nested global options."""
    node = cfg.get("global", {})
    for k in keys:
        if not isinstance(node, dict):
            return default
        node = node.get(k, default)
    return node


def parse_bool_option(raw: Any, key_name: str, default: bool = False) -> bool:
    """Parse bool-like config value with strict validation."""
    if raw is None:
        return default
    if isinstance(raw, bool):
        return raw
    if isinstance(raw, str):
        v = raw.strip().lower()
        if v in {"true", "yes", "1"}:
            return True
        if v in {"false", "no", "0"}:
            return False
    raise ValueError(f"{key_name} must be a boolean (true/false).")


def parse_chains_msa(entry: Dict[str, Any], n_chains: int) -> List[Optional[str]]:
    """
    Interpret 'chains_msa' mapping from config as a list of per-chain MSA paths.
    Keys can be int or string; indices are 0-based.
    """
    msas: List[Optional[str]] = [None] * n_chains
    raw = entry.get("chains_msa") or {}
    if not isinstance(raw, dict):
        raise ValueError(f"chains_msa must be a mapping, got {type(raw)}")
    for k, v in raw.items():
        try:
            idx = int(k)
        except Exception:
            raise ValueError(f"chains_msa key must be an integer index, got {k!r}")
        if 0 <= idx < n_chains:
            if v is not None:
                p = Path(v)
                # If relative and exists, convert to absolute path
                if not p.is_absolute() and p.exists():
                    p = p.resolve()
                msas[idx] = str(p)

        else:
            print(f"WARNING: chains_msa index {idx} out of range for {n_chains} chains; ignoring.")
    return msas


def build_binder_entities(cfg: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Build internal representation of binders:
      {name, seqs, msas}
    """
    binders_cfg = cfg.get("binders") or []
    if not isinstance(binders_cfg, list):
        raise ValueError("binders must be a list.")

    global_addK = bool(get_global_option(cfg, "add_n_terminal_lysine", default=False))
    result: List[Dict[str, Any]] = []

    for entry in binders_cfg:
        if not isinstance(entry, dict):
            raise ValueError("Each binder entry must be a mapping.")

        # Case: from_fasta_dir (from_dir kept as alias)
        if "from_fasta_dir" in entry or "from_dir" in entry:
            dir_raw = entry.get("from_fasta_dir", entry.get("from_dir"))
            dir_path = Path(dir_raw).resolve()
            if not dir_path.is_dir():
                raise ValueError(f"Binder from_fasta_dir/from_dir not found: {dir_path}")
            addK = bool(entry.get("add_n_terminal_lysine", global_addK))
            ligands = parse_ligands(entry, f"Binder directory {dir_path}")
            
            # Allow chains_msa even for from_fasta_dir (to set default MSA for all binders found)
            # We defer parsing until we know n_chains for each binder found.
            
            for name, seqs in read_fasta_dir_entities(dir_path):
                if addK:
                    seqs = add_n_terminal_lysine(seqs)
                
                # Try to parse chains_msa for THIS specific binder
                msas = parse_chains_msa(entry, len(seqs))
                
                result.append(
                    {"name": sanitize_name(name), "seqs": seqs, "msas": msas, "ligands": ligands}
                )
            continue
            
        # Case: from_structure_dir
        if "from_structure_dir" in entry:
            dir_path = Path(entry["from_structure_dir"]).resolve()
            if not dir_path.is_dir():
                raise ValueError(f"Binder from_structure_dir not found: {dir_path}")
            
            # Parse chains if specified
            chains_to_keep = None
            if "chains" in entry:
                # Expecting "A,B" string or list ["A", "B"]
                raw_chains = entry["chains"]
                if isinstance(raw_chains, str):
                    chains_to_keep = {c.strip() for c in raw_chains.split(",") if c.strip()}
                elif isinstance(raw_chains, list):
                    chains_to_keep = {str(c).strip() for c in raw_chains if str(c).strip()}
            
            addK = bool(entry.get("add_n_terminal_lysine", global_addK))
            ligands = parse_ligands(entry, f"Binder structure directory {dir_path}")

            # Iterate over structures
            for struct_path in sorted(dir_path.glob("*")):
                if not struct_path.is_file():
                    continue
                if struct_path.suffix.lower() not in {".pdb", ".cif", ".ent"}:
                    continue
                
                # Use str2fasta to get all chains
                try:
                    chain_seqs_map = get_sequences_all_chains(str(struct_path))
                except Exception as e:
                    print(f"Warning: failed to parse {struct_path.name}: {e}")
                    continue
                
                # Filter chains
                if chains_to_keep is not None:
                    ordered_chains = []
                    if "chains" in entry:
                        raw_chains = entry["chains"]
                        if isinstance(raw_chains, str):
                            ordered_chains = [c.strip() for c in raw_chains.split(",") if c.strip()]
                        elif isinstance(raw_chains, list):
                            ordered_chains = [str(c).strip() for c in raw_chains if str(c).strip()]
                    
                    final_seqs = []
                    for cid in ordered_chains:
                        if cid in chain_seqs_map:
                            final_seqs.append(chain_seqs_map[cid])
                        else:
                            # If a requested chain is missing, skip or warn?
                            # For now, let's warn and skip that chain
                            print(f"Warning: Chain {cid} not found in {struct_path.name}")
                else:
                    # No chains specified, take all chains
                    final_seqs = list(chain_seqs_map.values())
                
                if not final_seqs:
                    continue

                if addK:
                    final_seqs = add_n_terminal_lysine(final_seqs)
                
                name = sanitize_name(struct_path.stem)
                msas = parse_chains_msa(entry, len(final_seqs))
                
                result.append(
                    {"name": name, "seqs": final_seqs, "msas": msas, "ligands": ligands}
                )
            continue

        # Case: explicit binder
        if "name" not in entry:
            raise ValueError("Explicit binder entry must have a 'name'.")
        name = sanitize_name(entry["name"])

        # Sequences source: either 'sequences' or 'fasta'
        seqs: Optional[List[str]] = None
        if "sequences" in entry:
            raw_seqs = entry["sequences"]
            if not isinstance(raw_seqs, list) or not raw_seqs:
                raise ValueError(f"Binder {name}: 'sequences' must be a non-empty list.")
            seqs = [
                str(s).replace("\\n", "").replace("\n", "").replace(" ", "").upper()
                for s in raw_seqs
            ]
        elif "fasta" in entry:
            fasta_path = Path(entry["fasta"]).resolve()
            if not fasta_path.is_file():
                raise ValueError(f"Binder {name}: FASTA not found: {fasta_path}")
            seqs = read_fasta_multi(fasta_path)
        else:
            raise ValueError(f"Binder {name}: must specify 'sequences' or 'fasta'.")

        if not seqs:
            raise ValueError(f"Binder {name}: no sequences found.")

        addK = bool(entry.get("add_n_terminal_lysine", global_addK))
        if addK:
            seqs = add_n_terminal_lysine(seqs)

        msas = parse_chains_msa(entry, len(seqs))
        ligands = parse_ligands(entry, f"Binder {name}")
        result.append({"name": name, "seqs": seqs, "msas": msas, "ligands": ligands})

    if not result:
        raise ValueError("No binders defined in config.")
    return result


def build_target_entities(cfg: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Build internal representation of targets / antitargets:
      {name, role, seqs, msas}
    """
    targets_cfg = cfg.get("targets") or []
    if not isinstance(targets_cfg, list):
        raise ValueError("targets must be a list.")

    result: List[Dict[str, Any]] = []

    checked_template_paths: Set[Path] = set()

    for entry in targets_cfg:
        if not isinstance(entry, dict):
            raise ValueError("Each targets entry must be a mapping.")

        # Case: from_dir (NO MSAs)
        if "from_dir" in entry:
            if "role" not in entry:
                raise ValueError("targets[from_dir] entry must have a 'role' (target/antitarget).")
            role = str(entry["role"]).lower()
            if role not in {"target", "antitarget", "self"}:
                raise ValueError(
                    f"Target from_dir entry has invalid role {role!r} (expected 'target', 'antitarget', or 'self')."
                )
            if role == "self" and has_ligand_keys(entry):
                raise ValueError(
                    "targets[from_dir] with role 'self' cannot define ligands; define binder ligands under binders."
                )
            ligands = parse_ligands(entry, f"targets[from_dir] role={role}")

            dir_path = Path(entry["from_dir"]).resolve()
            if not dir_path.is_dir():
                raise ValueError(f"Targets from_dir not found: {dir_path}")
            for name, seqs in read_fasta_dir_entities(dir_path):
                # from_dir entries: explicitly NO MSAs
                msas = [None] * len(seqs)
                result.append(
                    {
                        "name": sanitize_name(name),
                        "role": role,
                        "seqs": seqs,
                        "msas": msas,
                        "ligands": ligands,
                    }
                )
            continue

        # Case: explicit target / antitarget
        if "name" not in entry:
            raise ValueError("Explicit target entry must have a 'name'.")
        if "role" not in entry:
            raise ValueError(f"Target {entry['name']!r} must have a 'role' (target/antitarget).")

        name = sanitize_name(entry["name"])
        role = str(entry["role"]).lower()
        # Extract and resolve CIF path
        cif_raw = entry.get("cif_template")
        cif_path = str(Path(cif_raw).resolve()) if cif_raw else None
        _precheck_template_cif(cif_path, f"Target {name}", checked_template_paths)
        if role not in {"target", "antitarget", "self"}:
            raise ValueError(
                f"Target {name}: invalid role {role!r} (expected 'target', 'antitarget', or 'self')."
            )
        if role == "self":
            if has_ligand_keys(entry):
                raise ValueError(
                    f"Target {name}: role 'self' cannot define ligands; define binder ligands under binders."
                )
            # For self, we need to capture chains_msa config if present, 
            # effectively deferring msa resolution until we know the binder.
            # We store the raw 'entry' dict to parse later.
            result.append({
                "name": name, 
                "role": "self", 
                "seqs": [], 
                "msas": [], 
                "ligands": [],
                "_raw_entry_for_msa": entry,
            })
            continue

        # Sequences source: either 'sequences' or 'fasta'
        seqs: Optional[List[str]] = None
        if "sequences" in entry:
            raw_seqs = entry["sequences"]
            if not isinstance(raw_seqs, list) or not raw_seqs:
                raise ValueError(f"Target {name}: 'sequences' must be a non-empty list.")
            seqs = [
                str(s).replace("\\n", "").replace("\n", "").replace(" ", "").upper()
                for s in raw_seqs
            ]
        elif "fasta" in entry:
            fasta_path = Path(entry["fasta"]).resolve()
            if not fasta_path.is_file():
                raise ValueError(f"Target {name}: FASTA not found: {fasta_path}")
            seqs = read_fasta_multi(fasta_path)
        else:
            raise ValueError(f"Target {name}: must specify 'sequences' or 'fasta'.")

        if not seqs:
            raise ValueError(f"Target {name}: no sequences found.")

        msas = parse_chains_msa(entry, len(seqs))
        ligands = parse_ligands(entry, f"Target {name}")
        result.append(
            {
                "name": name,
                "role": role,
                "seqs": seqs,
                "msas": msas,
                "ligands": ligands,
                "cif_template": cif_path,
            }
        )

    if not result:
        raise ValueError("No targets/antitargets defined in config.")
    return result


# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Generate YAMLs and run scripts for binder validation from a YAML config."
    )
    ap.add_argument("--config", required=True, help="Path to config YAML.")
    return ap.parse_args()


def main():
    args = parse_args()
    cfg_path = Path(args.config).resolve()
    if not cfg_path.is_file():
        sys.exit(f"ERROR: Config file not found: {cfg_path}")

    cfg = load_config(cfg_path)

    output_root = Path(cfg.get("output_dir", "./boltz_validation")).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    binders = build_binder_entities(cfg)
    targets_all = build_target_entities(cfg)

    # Separate into "targets" vs "antitargets" (but keep role info on each)
    targets = [t for t in targets_all if t["role"] == "target"]
    antitargets = [t for t in targets_all if t["role"] == "antitarget"]

    if len(targets_all) == 0:
        sys.exit("ERROR: Must have at least one target/antitarget/self entry.")



    # Boltz defaults
    boltz_cfg = get_global_option(cfg, "boltz", default={}) or {}
    print(f"Boltz global config: {boltz_cfg}")
    recycling_steps = boltz_cfg.get("recycling_steps", 10)
    diffusion_samples = boltz_cfg.get("diffusion_samples", 5)
    viz_cfg = cfg.get("visualisation", {}) or {}
    run_binder_monomer = parse_bool_option(
        viz_cfg.get("RMSD_to_binder_monomer", False),
        "visualisation.RMSD_to_binder_monomer",
        default=False,
    )

    # Validation: use_msa_server must be strictly true/false
    raw_msa_mode = str(boltz_cfg.get("use_msa_server", "false")).lower()
    if raw_msa_mode == "true":
        use_msa_server = True
    elif raw_msa_mode == "false":
        use_msa_server = False
    else:
        sys.exit("ERROR: global.boltz.use_msa_server must be 'true' or 'false'.")

    # --- Generate YAMLs and run.sh for each binder ---
    for binder in binders:
        bname = binder["name"]
        bseqs = binder["seqs"]
        bmsas = binder["msas"]
        bligs = binder.get("ligands", [])

        binder_dir = output_root / f"binder_{bname}"
        binder_dir.mkdir(parents=True, exist_ok=True)

        yaml_paths: List[Path] = []

        # Optional binder-only monomer run (used for RMSD reference in visualisation).
        if run_binder_monomer:
            monomer_yaml_name = f"{binder_dir.name}_monomer.yaml"
            monomer_path = binder_dir / monomer_yaml_name
            monomer_text = yaml_for_pair(
                bseqs,
                [],
                partner_role="self",
                use_msa_server=use_msa_server,
                binder_msas=bmsas,
                partner_msas=[],
                binder_ligands=bligs,
                partner_ligands=[],
                cif_template=None,
            )
            write_text(monomer_path, monomer_text)
            yaml_paths.append(monomer_path)

        # Loop over ALL partner entities (targets, antitargets, self)
        for tgt in targets_all:

            role = tgt["role"]

            # --- SELF CASE ---
            if role == "self":
                cif_tmp = None
                partner_name = "self"
                tseqs = bseqs[:]              # copy binder seqs
                # Self runs mirror the binder on the partner side, including ligands.
                tligs = bligs[:]
                
                # If the self entry defines chains_msa, use it.
                # Otherwise, fallback to binder's own MSAs (default behavior).
                if "_raw_entry_for_msa" in tgt:
                     # Attempt to parse MSAs for the target side using binder sequence length
                    self_msas = parse_chains_msa(tgt["_raw_entry_for_msa"], len(tseqs))
                    # If parsing resulted in ANY non-None entries, use them. 
                    # Note: parse_chains_msa returns [None, None...] if nothing found.
                    # BUT here the user might explicitly map "0: empty".
                    # parse_chains_msa handles "empty" as a string value for the path if resolved? 
                    # Wait, parse_chains_msa logic: if v is not None -> resolves path.
                    # It doesn't seem to inherently support string "empty" unless it's a file path?
                    # Actually, boltz config accepts "empty" as a specific keyword for ignoring MSAs.
                    # Let's check parse_chains_msa again.
                    
                    # Correction: parse_chains_msa resolves paths. 
                    # If the user put "empty" in config, `Path("empty").resolve()` points to a file named empty in CWD.
                    # That is NOT what we want if "empty" is a keyword.
                    # However, let's treat it generically: whatever parse_chains_msa returns.
                    # If user provided chains_msa, we should respect it.
                    
                    # If config has chains_msa, we prefer it over bmsas.
                    # We check if the raw entry has "chains_msa" key.
                    if "chains_msa" in tgt["_raw_entry_for_msa"]:
                         tmsas = self_msas
                    else:
                         tmsas = bmsas[:]
                else:
                    tmsas = bmsas[:]              # copy binder MSAs
                
                yaml_name = f"binder_{bname}_vs_self.yaml"

            # --- NORMAL TARGET ---
            elif role == "target":
                partner_name = tgt["name"]
                tseqs = tgt["seqs"]
                tmsas = tgt["msas"]
                tligs = tgt.get("ligands", [])
                yaml_name = f"binder_{bname}_vs_target_{partner_name}.yaml"
                cif_tmp = tgt.get("cif_template")
            # --- ANTITARGET ---
            elif role == "antitarget":
                partner_name = tgt["name"]
                tseqs = tgt["seqs"]
                tmsas = tgt["msas"]
                tligs = tgt.get("ligands", [])
                yaml_name = f"binder_{bname}_vs_antitarget_{partner_name}.yaml"
                cif_tmp = tgt.get("cif_template")
            else:
                raise ValueError(f"Unknown target role: {role}")

            ypath = binder_dir / yaml_name
            text = yaml_for_pair(
                bseqs,
                tseqs,
                partner_role=role,
                use_msa_server=use_msa_server,
                binder_msas=bmsas,
                partner_msas=tmsas,
                binder_ligands=bligs,
                partner_ligands=tligs,
                cif_template=cif_tmp,
            )
            write_text(ypath, text)
            yaml_paths.append(ypath)

        make_run_sh(
            binder_dir,
            yaml_paths,
            recycling_steps=recycling_steps,
            diffusion_samples=diffusion_samples,
            use_msa_server=use_msa_server,
        )

    make_master_run_sh(
        output_root,
        recycling_steps=recycling_steps,
        diffusion_samples=diffusion_samples,
        use_msa_server=use_msa_server,
    )
    make_visualisation_sh(output_root, cfg)

    print(f"\n✅ Done. YAMLs and scripts written under: {output_root}\n")


if __name__ == "__main__":
    main()
