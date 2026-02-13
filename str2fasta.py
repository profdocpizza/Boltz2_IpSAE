import os
from os import path as p
import argparse

from Bio.PDB import PDBParser, MMCIFParser, is_aa
from Bio.SeqUtils import seq1


def get_structure_parser(pdb_file):
    """Choose the correct parser based on file extension."""
    if pdb_file.endswith(".cif"):
        return MMCIFParser(QUIET=True)
    else:
        return PDBParser(QUIET=True)


def get_sequences_all_chains(structure_file):
    """
    Extract amino-acid sequences for all chains in a PDB or mmCIF file.
    Returns: dict {chain_id: sequence}
    """
    parser = get_structure_parser(structure_file)
    structure_id = p.basename(structure_file)
    structure = parser.get_structure(structure_id, structure_file)

    sequences = {}

    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                # Skip waters / heteroatoms
                if not is_aa(residue, standard=True):
                    continue

                resname = residue.get_resname()
                try:
                    aa = seq1(resname)
                    seq.append(aa)
                except Exception:
                    print(f"Warning: Unknown residue '{resname}' in {structure_file}")

            if seq:
                sequences[chain.id] = "".join(seq)

        # only use first model
        break

    return sequences


def main():
    args = parse_args()
    in_dir = args.indir
    chain_id = args.chain

    structure_files = [
        p.join(in_dir, f)
        for f in os.listdir(in_dir)
        if f.endswith(".pdb") or f.endswith(".cif")
    ]

    for structure_file in structure_files:
        base_name = p.basename(structure_file).rsplit(".", 1)[0]
        fasta_name = p.join(in_dir, base_name + ".fasta")

        sequences = get_sequences_all_chains(structure_file)

        with open(fasta_name, "w") as f:
            if chain_id is None:
                for cid, seq in sequences.items():
                    f.write(f">{base_name}_{cid}\n{seq}\n")
            else:
                if chain_id in sequences:
                    f.write(f">{base_name}\n{sequences[chain_id]}\n")
                else:
                    print(f"⚠ Chain {chain_id} not found in {structure_file}")

        print(f"✅ Saved {fasta_name}")

    if args.merge:
        merged_path = p.join(in_dir, "combined.fasta")
        fasta_files = [
            p.join(in_dir, f) for f in os.listdir(in_dir) if f.endswith(".fasta")
        ]

        with open(merged_path, "w") as outfile:
            for ff in fasta_files:
                with open(ff) as infile:
                    outfile.write(infile.read().strip() + "\n")

        print(f"📘 Merged {len(fasta_files)} FASTA files → {merged_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences from PDB or mmCIF files in a directory."
    )
    parser.add_argument(
        "--indir", "-i", required=True,
        help="Directory containing .pdb or .cif files."
    )
    parser.add_argument(
        "--chain", "-c", default=None,
        help="Chain ID to extract. Default: extract ALL chains."
    )
    parser.add_argument(
        "--merge", "-m", action="store_true",
        help="Merge all generated FASTA files into a single multi-FASTA file."
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
