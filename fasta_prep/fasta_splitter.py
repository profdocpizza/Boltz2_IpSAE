#!/usr/bin/env python3
"""
fasta_splitter.py
Split a multi-FASTA file into individual FASTA files based on entry names.

Usage:
    python fasta_splitter.py input.fasta output_dir/
"""

import os
import sys
from pathlib import Path

def split_fasta(input_fasta, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_fasta, 'r') as f:
        current_header = None
        current_seq_lines = []

        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                # write previous entry
                if current_header:
                    write_entry(current_header, current_seq_lines, output_dir)
                current_header = line[1:].split()[0]  # take first token as filename
                current_seq_lines = []
            else:
                current_seq_lines.append(line)

        # write last entry
        if current_header:
            write_entry(current_header, current_seq_lines, output_dir)


def write_entry(header, seq_lines, output_dir):
    filename = f"{header}.fasta"
    filepath = output_dir / filename
    with open(filepath, 'w') as out:
        out.write(f">{header}\n")
        out.write("\n".join(seq_lines) + "\n")


def main():
    if len(sys.argv) != 3:
        print("Usage: python fasta_splitter.py input.fasta output_dir/")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    split_fasta(input_fasta, output_dir)


if __name__ == "__main__":
    main()
