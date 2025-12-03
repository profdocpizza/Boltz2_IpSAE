#!/usr/bin/env python3
import os
import sys
from pathlib import Path

def merge_fasta(directory):
    directory = Path(directory)
    if not directory.is_dir():
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    fasta_files = sorted(
        [f for f in directory.iterdir() if f.suffix.lower() in [".fasta", ".fa"]]
    )

    if not fasta_files:
        print("No FASTA files found.")
        sys.exit(0)

    output_path = directory / "merged.fasta"

    with open(output_path, "w") as outfile:
        for fasta in fasta_files:
            with open(fasta, "r") as infile:
                outfile.write(infile.read().rstrip() + "\n")

    print(f"Merged {len(fasta_files)} FASTA files into: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <directory>")
        sys.exit(1)

    merge_fasta(sys.argv[1])
