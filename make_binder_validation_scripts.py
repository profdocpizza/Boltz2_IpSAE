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
  - MSA can be provided for any chain via config (chains_msa).
  - from_dir entries NEVER have MSAs (by design).
  - Uses 'target_' / 'antitarget_' prefixes in YAML names:
       binder_<binder>_vs_target_<name>.yaml
       binder_<binder>_vs_antitarget_<name>.yaml
  - Generates:
       - Per-binder YAMLs for all binder–(anti)target pairs
       - Per-binder run.sh
       - Global run_all_cofolding.sh
       - Visualization helper script
"""

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

try:
    import yaml
except ImportError:
    sys.exit("ERROR: PyYAML is required. Install with `pip install pyyaml`.")

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



def yaml_for_pair(
    binder_seqs: List[str],
    partner_seqs: List[str],
    partner_role: str,
    binder_msas: Optional[List[Optional[str]]] = None,
    partner_msas: Optional[List[Optional[str]]] = None,
) -> str:
    """
    Build Boltz YAML for a binder–partner pair.

    - Binder chains first with IDs: A, B, C, ...
    - Partner chains next with IDs depending on role:
          target     → TA, TB, ...
          antitarget → AA, AB, ...
          other      → PA, PB, ...

    If any MSA is provided for binder or partner chains, then EVERY chain gets
    an 'msa:' entry. Chains without a file get 'msa: empty'.

    YAML format for Boltz stays:

        version: 1
        sequences:
          - protein:
              id: ...
              sequence: ...
              msa: <optional field>
    """

    lines: List[str] = ["version: 1", "sequences:"]

    binder_msas = binder_msas or [None] * len(binder_seqs)
    partner_msas = partner_msas or [None] * len(partner_seqs)

    # --- Binder chains (A, B, ...) ---
    for i, seq in enumerate(binder_seqs):
        cid = chr(ord("A") + i)
        lines.append("  - protein:")
        lines.append(f"      id: {cid}")
        lines.append(f"      sequence: {seq}")
        msa_path = binder_msas[i] if i < len(binder_msas) and binder_msas[i] else "empty"
        lines.append(f"      msa: {msa_path}")

    # --- Partner chains (TA/TB/... or AA/AB/...) ---
    for i, seq in enumerate(partner_seqs):
        cid = _partner_chain_id(partner_role, i)
        lines.append("  - protein:")
        lines.append(f"      id: {cid}")
        lines.append(f"      sequence: {seq}")
        msa_path = partner_msas[i] if i < len(partner_msas) and partner_msas[i] else "empty"
        lines.append(f"      msa: {msa_path}")

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


def make_master_run_sh(output_root: Path) -> None:
    """Generate top-level script to run all binder run.sh scripts with paths relative to the script itself."""
    lines = [
        "#!/bin/bash",
        "set -e",
        "",
        '# Determine directory of this script',
        'DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"',
        "",
        "# Run all binder_* run.sh scripts relative to script location",
        'for f in $(find "$DIR" -type f -name "run.sh" | sort); do',
        '  echo "Running $f..."',
        '  (cd "$(dirname "$f")" && bash run.sh)',
        "done",
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

    # Flags
    best_flag = "--use_best_model" if use_best else ""
    cpu_flag = f"--num_cpu {num_cpu}" if num_cpu is not None else ""

    script_line = (
        f"python {os.path.dirname(os.path.abspath(__file__))}/visualise_binder_validation.py "
        f"--ipsae_e {ipsae_e} "
        f"--ipsae_d {ipsae_d} "
        f"{cpu_flag} "
        f"--root_dir {output_root} "
        f"--generate_data --plot {best_flag}"
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
                msas[idx] = str(v)
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

        # Case: from_dir (NO MSAs)
        if "from_dir" in entry:
            dir_path = Path(entry["from_dir"]).resolve()
            if not dir_path.is_dir():
                raise ValueError(f"Binder from_dir not found: {dir_path}")
            addK = bool(entry.get("add_n_terminal_lysine", global_addK))
            for name, seqs in read_fasta_dir_entities(dir_path):
                if addK:
                    seqs = add_n_terminal_lysine(seqs)
                # from_dir entries: explicitly NO MSAs
                msas = [None] * len(seqs)
                result.append(
                    {"name": sanitize_name(name), "seqs": seqs, "msas": msas}
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
        result.append({"name": name, "seqs": seqs, "msas": msas})

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
                    f"Target {name}: invalid role {role!r} (expected 'target', 'antitarget', or 'self')."
                )

            dir_path = Path(entry["from_dir"]).resolve()
            if not dir_path.is_dir():
                raise ValueError(f"Targets from_dir not found: {dir_path}")
            for name, seqs in read_fasta_dir_entities(dir_path):
                # from_dir entries: explicitly NO MSAs
                msas = [None] * len(seqs)
                result.append(
                    {"name": sanitize_name(name), "role": role, "seqs": seqs, "msas": msas}
                )
            continue

        # Case: explicit target / antitarget
        if "name" not in entry:
            raise ValueError("Explicit target entry must have a 'name'.")
        if "role" not in entry:
            raise ValueError(f"Target {entry['name']!r} must have a 'role' (target/antitarget).")

        name = sanitize_name(entry["name"])
        role = str(entry["role"]).lower()
        if role not in {"target", "antitarget", "self"}:
            raise ValueError(
                f"Target {name}: invalid role {role!r} (expected 'target', 'antitarget', or 'self')."
            )
        if role == "self":
            seqs = []   # placeholder to indicate "binder will fill this"
            msas = []   # placeholder for binder MSAs
            result.append({"name": name, "role": "self", "seqs": seqs, "msas": msas})
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
        result.append({"name": name, "role": role, "seqs": seqs, "msas": msas})

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
    use_msa_server_mode = str(boltz_cfg.get("use_msa_server", "auto")).lower()
    if use_msa_server_mode not in {"auto", "true", "false"}:
        sys.exit("ERROR: global.boltz.use_msa_server must be 'auto', 'true', or 'false'.")


    # --- Generate YAMLs and run.sh for each binder ---
    for binder in binders:
        bname = binder["name"]
        bseqs = binder["seqs"]
        bmsas = binder["msas"]

        binder_dir = output_root / f"binder_{bname}"
        binder_dir.mkdir(parents=True, exist_ok=True)

        yaml_paths: List[Path] = []
        any_pair_uses_msa = False

        # Loop over ALL partner entities (targets, antitargets, self)
        for tgt in targets_all:

            role = tgt["role"]

            # --- SELF CASE ---
            if role == "self":
                partner_name = "self"
                tseqs = bseqs[:]              # copy binder seqs
                tmsas = bmsas[:]              # copy binder MSAs
                yaml_name = f"binder_{bname}_vs_self.yaml"

            # --- NORMAL TARGET ---
            elif role == "target":
                partner_name = tgt["name"]
                tseqs = tgt["seqs"]
                tmsas = tgt["msas"]
                yaml_name = f"binder_{bname}_vs_target_{partner_name}.yaml"

            # --- ANTITARGET ---
            elif role == "antitarget":
                partner_name = tgt["name"]
                tseqs = tgt["seqs"]
                tmsas = tgt["msas"]
                yaml_name = f"binder_{bname}_vs_antitarget_{partner_name}.yaml"

            else:
                raise ValueError(f"Unknown target role: {role}")

            ypath = binder_dir / yaml_name
            text = yaml_for_pair(
                bseqs,
                tseqs,
                partner_role=role,
                binder_msas=bmsas,
                partner_msas=tmsas,
            )
            write_text(ypath, text)
            yaml_paths.append(ypath)

            if any(bmsas) or any(tmsas):
                any_pair_uses_msa = True

        # Decide whether to use MSA server
        if use_msa_server_mode == "true":
            use_msa_server = True
        elif use_msa_server_mode == "false":
            use_msa_server = False
        else:
            use_msa_server = not any_pair_uses_msa

        make_run_sh(
            binder_dir,
            yaml_paths,
            recycling_steps=recycling_steps,
            diffusion_samples=diffusion_samples,
            use_msa_server=use_msa_server,
        )

    make_master_run_sh(output_root)
    make_visualisation_sh(output_root, cfg)

    print(f"\n✅ Done. YAMLs and scripts written under: {output_root}\n")


if __name__ == "__main__":
    main()
