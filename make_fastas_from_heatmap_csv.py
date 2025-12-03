import argparse
import os
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Filter binder sequences by thresholds.")
    parser.add_argument("--csv_file", required=True,
                        help="Path to ipSAE_min_heatmap.csv")
    parser.add_argument("--output_dir", default=".",
                        help="Directory to write FASTA files")
    parser.add_argument("--self_threshold", type=float, default=0.15)
    parser.add_argument("--target_threshold", type=float, default=0.3)
    parser.add_argument("--antitarget_threshold", type=float, default=0.15)

    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)

    # Apply thresholds
    mask = (
        (df["self"] < args.self_threshold) &
        (df["target"] > args.target_threshold) &
        (df["antitarget"] < args.antitarget_threshold)
    )

    filtered = df[mask]

    if filtered.empty:
        print("No binders passed thresholds.")
        return

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Write each binder as its own FASTA file
    for _, row in filtered.iterrows():
        binder_name = str(row["binder_short"]).strip()
        sequence = row["binder_sequence"]

        fasta_path = os.path.join(args.output_dir, f"{binder_name}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{binder_name}\n{sequence}\n")

        print(f"Saved: {fasta_path}")

    print(f"\nDone! {len(filtered)} FASTA files written.")

if __name__ == "__main__":
    main()
