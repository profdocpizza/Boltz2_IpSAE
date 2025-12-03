import os
from os import path as p
from pdbUtils.pdbUtils import pdb2df
import argparse


def init_3_to_1():
    return {
        "ALA": "A","ARG": "R","ASN": "N","ASP": "D","CYS": "C",
        "GLN": "Q","GLU": "E","GLY": "G","HIS": "H","ILE": "I",
        "LEU": "L","LYS": "K","MET": "M","PHE": "F","PRO": "P",
        "SER": "S","THR": "T","TRP": "W","TYR": "Y","VAL": "V"
    }


def get_sequences_all_chains(pdb_file):
    three_to_one = init_3_to_1()
    df = pdb2df(pdb_file)
    sequences = {}

    for cid, chain_df in df.groupby("CHAIN_ID"):
        seq = ""
        for _, res_df in chain_df.groupby("RES_ID"):
            res_name = res_df.iloc[0]["RES_NAME"]
            aa = three_to_one.get(res_name)
            if aa:
                seq += aa
            else:
                print(f"Warning: Unknown residue '{res_name}' in {pdb_file}")
        sequences[cid] = seq

    return sequences


def main():
    args = parse_args()
    in_dir = args.indir
    chain_id = args.chain

    pdb_files = [p.join(in_dir, f) for f in os.listdir(in_dir) if f.endswith(".pdb")]

    for pdb_file in pdb_files:
        pdb_name = p.basename(pdb_file).split(".")[0]
        fasta_name = p.join(in_dir, pdb_name + ".fasta")

        with open(fasta_name, "w") as f:
            sequences = get_sequences_all_chains(pdb_file)

            if chain_id is None:
                for cid, seq in sequences.items():
                    f.write(f">{pdb_name}_{cid}\n{seq}\n")
            else:
                if chain_id in sequences:
                    f.write(f">{pdb_name}\n{sequences[chain_id]}\n")
                else:
                    print(f"⚠ Chain {chain_id} not found in {pdb_file}")

        print(f"✅ Saved {fasta_name}")

    if args.merge:
        merged_path = p.join(in_dir, "combined.fasta")
        fasta_files = [p.join(in_dir, f) for f in os.listdir(in_dir) if f.endswith(".fasta")]

        with open(merged_path, "w") as outfile:
            for ff in fasta_files:
                with open(ff) as infile:
                    outfile.write(infile.read().strip() + "\n")

        print(f"📘 Merged {len(fasta_files)} FASTA files → {merged_path}")



def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences from PDB files in a directory."
    )
    parser.add_argument(
        "--indir", "-i", required=True, help="Directory containing .pdb files."
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
