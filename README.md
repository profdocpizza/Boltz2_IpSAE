# Boltz2_IpSAE

This repository contains scripts for co-folding using Boltz2, calculating and visualising [ipSAE](https://github.com/DunbrackLab/IPSAE) metric. With a single .yaml file you can cofold your binder to multiple targets and itself to monitor specificity of binding.

![Example output](images/ipSAE_min_heatmap.png)


## Configuration

The `config.yaml` file controls all aspects of binder validation. Below is a comprehensive reference with all supported fields.

### Config Structure Overview

```yaml
version: 1                                              # Required. Must be 1.
output_dir: ./boltz_validation                         # Required. Where all outputs will be written.

global:                                                 # Optional. Global settings for all predictions.
  boltz:                                                # Optional. Boltz2 parameters.
    recycling_steps: 10                                 # Default: 10. Number of recycling steps in Boltz.
    diffusion_samples: 5                                # Default: 5. Number of diffusion samples.
    use_msa_server: true                                # Default: false. Use ColabFold MSA server. MUST be "true"/"false" (string).
  add_n_terminal_lysine: false                          # Default: false. Prepend 'K' to all binder sequences if missing.

visualisation:                                          # Optional. Parameters for visualization and analysis.
  ipsae_error_threshold: 15                             # Default: 15. ipSAE error threshold for plots.
  ipsae_distance_threshold: 15                          # Default: 15. ipSAE distance threshold for plots.
  use_best_model: false                                 # Default: true. If false, plots average across all models.
  num_cpu: 50                                           # Optional. Number of CPUs for parallel visualization. Omit to auto-detect.
  RMSD_to_binder_monomer: false                         # Default: false. Generate monomer predictions and compute RMSD.

binders:                                                # Required. List of binders to validate.
  - name: BinderA                                       # Optional. Name of binder. If using from_fasta_dir, auto-derived.
    sequences:                                          # Optional. Inline protein sequences (mutually exclusive with fasta, from_fasta_dir, from_structure_dir).
      - MKVLSPADKTNV...                                 # Each item is a single-chain or multi-line FASTA sequence.
    # OR use one of these alternatives:
    # fasta: /path/to/binder.fasta                      # Single FASTA file with one or more sequences (protein sequences in this file become chains A, B, C...).
    # from_fasta_dir: /path/to/binder_fasta_dir/        # Directory with multiple .fasta/.fa files (each file becomes one binder).
    # from_structure_dir: /path/to/structures/          # Directory with PDB/CIF files. Requires additional options below if used.
    #   chain_extraction:
    #     extract_method: select_chains
    #     chain_ids: [A, B]                             # Which chains to extract from each structure.

    chains_msa:                                         # Optional. Provide multiple sequence alignments per chain.
      # 0: empty                                        # Chain 0 (A): "empty" means skip MSA for this chain.
      # 1: /path/to/alignment.a3m                       # Chain 1 (B): path to MSA file.
    
    # Ligands for this binder (applied to all target/antitarget/self runs):
    # ligand_ccd:                                        # Optional. PDB chemical component codes (e.g., ZN for zinc).
    #   - ZN                                             # Single code or list of codes. Repeated codes become multiple chain IDs.
    #   - ATP
    # ligand_smiles:                                     # Optional. SMILES strings for small molecules.
    #   - "CCO"                                          # Single SMILES or list of SMILES strings.
    #   - "c1ccccc1"

targets:                                                # Required. List of target/antitarget partners (or self).
  - name: Target1                                       # Required. Name of target.
    role: target                                        # Required. One of: "target" (desired binding), "antitarget" (avoid binding), "self" (binder monomer).
    sequences:                                          # Optional. Inline sequences (mutually exclusive with fasta, from_dir).
      - EVQLQQSGPVLVK...                                # Heavy chain example.
      - DVLMTQTPLSLPV...                                # Light chain example.
    # OR use one of these alternatives:
    # fasta: /path/to/target.fasta                      # Single FASTA file with one or more sequences.
    # from_dir: /path/to/targets/                       # Directory with multiple .fasta files (requires role, only for targets).

    chains_msa:                                         # Optional. MSA per chain (same structure as binders).
      # 0: empty
      # 1: /path/to/heavy_chain.a3m
    
    cif_template: /path/to/template.cif                 # Optional. CIF structure file to guide co-folding (usually a known complex).
                                                        # If set, Boltz YAML will include this template without chain_id (Boltz auto-matches).
    
    # Ligands for this target (only for role: "target" or "antitarget", NOT for "self"):
    # ligand_ccd:                                        # Optional. Same as binder ligands.
    #   - FAD
    # ligand_smiles:                                     # Optional. Same as binder ligands.
    #   - "C1=CC=CC=C1"

  - name: self                                          # Special case: binder-only monomer prediction.
    role: self                                          # Required. Role "self" means binder predicts itself (no partner).
    chains_msa:                                         # Optional. MSA for binder chains (if different from binder's default).
      # 0: empty
```

### Key Rules and Constraints

**Binders:**
- Must define sequences via exactly ONE method: `sequences`, `fasta`, `from_fasta_dir`, or `from_structure_dir`
- `sequences` must be a list of strings (one per chain)
- `fasta` is a path to a single FASTA file; multiple sequences in the file become chains A, B, C...
- `from_fasta_dir` reads all `.fasta`/`.fa` files in a directory; each file name (stem) becomes a binder name
- `name` is required if using `sequences` or `fasta`; it is auto-derived from filename if using `from_fasta_dir`
- `chains_msa` is a dict mapping chain index (0, 1, 2...) to either a file path or "empty"

**Targets:**
- Must define sequences via exactly ONE method: `sequences`, `fasta`, or `from_dir`
- `role` is required and must be one of: `target`, `antitarget`, `self`
- `role: self` mirrors the binder against itself; it ignores any `sequences`, `fasta`, `from_dir` values (uses binder's sequences on both sides)
- `cif_template` is optional and recommended for known complexes; it helps guide structure prediction
- Ligands (`ligand_ccd`, `ligand_smiles`) are allowed ONLY for `role: target` or `role: antitarget`
  - `role: self` cannot define ligands; self runs duplicate the binder ligands onto both mirrored sides

**Global Options:**
- `use_msa_server` must be the string `"true"` or `"false"` (not a boolean)
- `add_n_terminal_lysine: true` will prepend 'K' to every binder sequence that doesn't already start with 'K'

**Visualization:**
- `RMSD_to_binder_monomer: true` generates an extra binder-only monomer run per binder and computes RMSD to that monomer

### Minimal Working Example

```yaml
version: 1
output_dir: ./boltz_validation

binders:
  - name: MyBinder
    sequences:
      - MKVLSPADKTNV

targets:
  - name: MyTarget
    role: target
    sequences:
      - EVQLQQSGPVLVK
```

### Practical Examples

**Example 1: Multiple binders from FASTA files**
```yaml
version: 1
output_dir: ./outputs

binders:
  - from_fasta_dir: ./binder_sequences/    # Each .fasta file becomes one binder (WT, M100K, M100A...).

targets:
  - name: Target9D9
    role: target
    sequences:
      - EVQLQQSGPVLVKPGASVKMSCKASGYTFTDYYMNWVK...
    cif_template: /path/to/9d9_template.cif
```

**Example 2: Binder with MSA and ligands**
```yaml
version: 1
output_dir: ./outputs

binders:
  - name: ProteinWithLigand
    sequences:
      - MKVLSPADKTNVVLWAGSKQ
    chains_msa:
      0: /path/to/alignment.a3m
    ligand_ccd:
      - ZN
      - ATP

targets:
  - name: TargetA
    role: target
    sequences:
      - EVQLQQSGPVLVKPGA
    ligand_ccd:
      - FAD
  
  - name: TargetB
    role: antitarget
    sequences:
      - MGLAILIFVTVLLISDAVSVETQAY
  
  - name: self
    role: self
```

**Example 3: Multi-chain target with CIF template**
```yaml
version: 1
output_dir: ./outputs

binders:
  - name: MyBinder
    sequences:
      - MKVLSPADKTNVVLWAGSKQ

targets:
  - name: ComplexTarget
    role: target
    sequences:
      - EVQLQQSGPVLVKPGA          # Heavy chain
      - DVLMTQTPLSLPVSLGDQ        # Light chain
    chains_msa:
      0: /path/to/hc.a3m
      1: /path/to/lc.a3m
    cif_template: /path/to/complex.cif
```

### What is NOT Supported

- You **cannot** have both `sequences` and `fasta` in the same binder/target entry
- `from_dir` for binders is **not supported** (use `from_fasta_dir` instead)
- MSA files for `role: self` targets are ignored; self always uses the binder's own chains_msa
- Ligands on `role: self` entries are **not allowed** (define them on the binder instead)
- Chain IDs in `cif_template` are **not supported**; Boltz auto-matches template chains to sequences


## How to Run

1.  **Clone the repository and set up the [`Boltz2`](https://github.com/jwohlwend/boltz) environment.**
2.  **Navigate to the example directory:**

    ```bash
    cd ./Boltz2_IpSAE/example_run
    ```

3.  **Run the main script to generate the validation scripts:**

    ```bash
    python ../make_binder_validation_scripts.py --config config.yaml
    ```

4.  **Run all co-folding validation scripts with a single script**

    ```bash
    ./boltz_validation/run_all_cofolding.sh
    ```

5.  **When co-folding is done, visualize the results:**

    ```bash
    ./boltz_validation/visualise_cofolding_results.sh
    ```

This will generate a plot of the results, which will be saved in the `boltz_validation/summary` directory.

## Ligands in config

Ligands can be added to binders and targets via the config file (see **Configuration** section above for details and examples).

**For binders:**
- `binders[*].ligand_ccd` - PDB chemical component codes (single string or list)
- `binders[*].ligand_smiles` - SMILES strings (single string or list)
- Ligands defined here are applied to all target/antitarget/self runs with this binder
  - For `role: self`, the binder ligands are duplicated so both mirrored copies of the binder carry them

**For targets:**
- `targets[*].ligand_ccd` / `targets[*].ligand_smiles` - **only allowed for `role: target` or `role: antitarget`**
- **Not allowed for `role: self`** - self runs duplicate ligands from the binder entry instead

**Behavior:**
- Each key accepts a single string or a list of strings
- Repeated ligand values are automatically grouped into one Boltz ligand entity with multiple chain IDs (e.g., `id: [LA, LB]`)
- For example: `ligand_ccd: [ZN, ZN]` creates one zinc entity with chains LA and LB

## Template CIF behavior

The `cif_template` field (optional, for targets) guides Boltz structure prediction with known complex structures:

- **When set:** `targets[*].cif_template: /path/to/template.cif`
  - The generated Boltz YAML includes only:
    ```yaml
    templates:
      - cif: /path/to/template.cif
    ```
  - **No `chain_id` is written** intentionally, so Boltz auto-selects and matches template chains to the predicted protein sequences
  
- **When not set:** Boltz predicts structure de novo from MSA and language model information

- **Precheck during config parsing:**
  - The script performs an early validation check on the CIF file (using Boltz parser when available, with Gemmi fallback)
  - If the template cannot be parsed, a warning is printed with a link to the PDB CIF converter:
    - https://mmcif.pdbj.org/converter/index.php?l=en
  
- **Common use case:** Provide the PDB structure of a known binder-target complex to improve co-folding accuracy

## Optional monomer RMSD column

Set `visualisation.RMSD_to_binder_monomer: true` to:

- generate an extra binder-only monomer Boltz prediction per binder (same global Boltz params as pair runs),
- compute RMSD from each binder-vs-* model to all binder monomer models (average across monomer samples),
- append this value to existing rows as `RMSD_to_binder_monomer` in `summary/ipsae_summary_all_binders.csv`.
