#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Visualize Boltz/ipSAE validation outputs.

Key assumptions / behaviour:

- Boltz configs use chain IDs:
    * Binder chains: A, B, C, ...
    * Target chains: TA, TB, TC, ...
    * Antitarget chains: AA, AB, AC, ...
- Binder is always the A-chain (chain_of_focus = "A").
- Boltz outputs live under:
      binder_<binder_name>/outputs/boltz_results_<yaml_stem>/
- YAML stems look like:
      binder_<binder>_vs_target_<target>
      binder_<binder>_vs_antitarget_<name>

This script:
  * runs ipSAE on all models for each binder–(anti)target pair
  * extracts metrics for chain A vs its best partner
  * stores:
        - binder
        - vs (full name)
        - partner (e.g. Spike, HA, ...)
        - target_type (target / antitarget / unknown)
        - model_idx
        - numeric ipSAE metrics (_min, _max)
  * makes per-binder stripplots (all targets & antitargets together)
  * makes global heatmaps (ipSAE_min, ipSAE_max) averaged across models
"""

import argparse
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path


def classify(chain_id):
    # Binder chains: A, B, C ...
    if len(chain_id) == 1 and chain_id.isupper():
        return "binder"

    # Self-chains: SA, SB, SC ...
    if chain_id.startswith("S"):
        return "self"

    # Targets: TA, TB, ...
    if chain_id.startswith("T"):
        return "target"

    # Antitarget: AA, AB, ...
    if chain_id.startswith("A") and len(chain_id) > 1:
        return "antitarget"

    return "other"



def run_ipsae(
    pae_file,
    cif_file,
    pae_cutoff=15,
    dist_cutoff=15,
    chain_of_focus="A"
):
    """
    Run ipsae.py and compute min/max metrics for chain_of_focus.

    If multiple partner chains exist, choose the one with the highest ipSAE_max
    (using only ASYM rows). Then compute metric_min / metric_max from only
    those rows.
    """

    import pandas as pd
    import numpy as np
    import os
    import subprocess

    # ---------------------------
    # Run ipSAE
    # ---------------------------
    cmd = [
        "python", f"{os.path.dirname(os.path.abspath(__file__))}/ipsae.py",
        str(pae_file),
        str(cif_file),
        str(pae_cutoff),
        str(dist_cutoff)
    ]
    subprocess.run(cmd, check=True)

    out_txt = str(cif_file).replace(".cif", f"_{pae_cutoff}_{dist_cutoff}.txt")
    if not os.path.exists(out_txt):
        raise FileNotFoundError(f"Missing ipSAE output: {out_txt}.\nCommand was: {' '.join(cmd)}")

    # ---------------------------
    # Load table

    # ---------------------------------------------------
    # Convert fixed-width text to CSV by collapsing spaces
    # ---------------------------------------------------
    with open(out_txt, "r") as f:
        raw_lines = f.readlines()

    clean_lines = []
    for line in raw_lines:
        stripped = line.strip()

        # Skip fully blank lines
        if not stripped:
            continue

        # Replace 2+ spaces with a single comma
        cleaned = re.sub(r"\s+", ",", stripped)

        clean_lines.append(cleaned)

    csv_tmp = out_txt + ".csv"

    # Write cleaned CSV
    with open(csv_tmp, "w") as f:
        for cl in clean_lines:
            f.write(cl + "\n")

    # Load CSV normally
    df = pd.read_csv(csv_tmp)
        

    df["Type"] = df["Type"].astype(str).str.lower().str.strip()
    df = df[df["Type"].str.contains("asym", na=False)]

    # Required columns
    chn1_col = "Chn1"
    chn2_col = "Chn2"

    # Compute category for each chain
    df["cat1"] = df[chn1_col].apply(classify)
    df["cat2"] = df[chn2_col].apply(classify)

    # Keep only binder–target or binder–antitarget pairs
    df = df[
        ((df["cat1"] == "binder") & (df["cat2"].isin(["target", "antitarget","self"]))) |
        ((df["cat2"] == "binder") & (df["cat1"].isin(["target", "antitarget","self"])))
    ]

    if df.empty:
        raise ValueError(f"No valid binder–target/antitarget ASYM rows for chain {chain_of_focus}")



    # ---------------------------
    # Keep only rows where focus chain appears
    # ---------------------------
    df = df[(df[chn1_col] == chain_of_focus) | (df[chn2_col] == chain_of_focus)]
    if df.empty:
        raise ValueError(f"No ASYM rows involving chain {chain_of_focus}")

    # ---------------------------
    # Identify partner chains
    # ---------------------------
    partners = set()
    for _, row in df.iterrows():
        partner = row[chn2_col] if row[chn1_col] == chain_of_focus else row[chn1_col]
        partners.add(partner)
    partners = sorted(partners)

    # ---------------------------
    # If multiple partners → choose highest ipSAE_max partner
    # ---------------------------
    partner_best = None
    partner_best_score = -np.inf
    best_df = None

    for p in partners:
        sub = df[
            ((df[chn1_col] == chain_of_focus) & (df[chn2_col] == p)) |
            ((df[chn2_col] == chain_of_focus) & (df[chn1_col] == p))
        ]
        if sub.empty:
            continue

        ipSAE_max = sub["ipSAE"].max()
        if ipSAE_max > partner_best_score:
            partner_best_score = ipSAE_max
            partner_best = p
            best_df = sub.copy()

    if partner_best is None or best_df is None:
        raise ValueError(f"No valid partner rows found for {chain_of_focus}")

    # ---------------------------
    # Compute min/max metrics from best partner rows only
    # ---------------------------
    numeric_cols = best_df.select_dtypes(include=[np.number]).columns.tolist()
    numeric_cols = [c for c in numeric_cols if c.lower() != "model"]

    output = {
        "chain_of_focus": chain_of_focus,
        "involved_chains": partner_best
    }

    for col in numeric_cols:
        output[f"{col}_min"] = best_df[col].min()
        output[f"{col}_max"] = best_df[col].max()

    return output


def parse_vs_name(vs_name: str):
    """
    Parse:
        binder_<binder>_vs_target_<name>
        binder_<binder>_vs_antitarget_<name>
        binder_<binder>_vs_self
    """
    m = re.search(r"_vs_(target|antitarget|self)_(.*)$", vs_name)
    if m:
        role = m.group(1)
        partner = m.group(2)
        return partner, role

    # self without a partner name (binder_X_vs_self)
    m2 = re.search(r"_vs_(self)$", vs_name)
    if m2:
        return "self", "self"

    return vs_name, "unknown"



def analyse_binder(binder_dir: Path ,args):
    """
    Analyse a binder directory: compute ipSAE for all vs_* pairs, save plots.
    """
    plots_dir = binder_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    binder_records = []

    for vs_dir in (binder_dir / "outputs").glob("boltz_results_*vs*"):
        vs_name = vs_dir.name.replace("boltz_results_", "")
        pred_root = vs_dir / "predictions" / vs_name

        partner_name, target_type = parse_vs_name(vs_name)

        for pae_file in pred_root.glob("pae_*_model_*.npz"):
            m = re.search(r"model_(\d+)", pae_file.name)
            if not m:
                continue
            model_idx = m.group(1)
            cif_file = pred_root / pae_file.name.replace("pae_", "").replace(".npz", ".cif")
            if not cif_file.exists():
                continue

            try:
                rec = run_ipsae(pae_file, cif_file, chain_of_focus="A",pae_cutoff=int(args.ipsae_e), dist_cutoff=int(args.ipsae_d))
            except Exception as e:
                print(f"⚠️ ipSAE failed for {pae_file} ({e}). Skipping.")
                continue
            # -----------------------------------------------------------
            # Add sequences from YAML (binder1, binder2, target1, ...)
            # -----------------------------------------------------------
            yaml_path = binder_dir / f"{vs_name}.yaml"

            binder_seqs, target_seqs, antitarget_seqs, self_seqs = extract_sequences_from_yaml(yaml_path)

            # binder may have multiple chains
            if len(binder_seqs) <= 1:
                rec["binder_sequence"] = binder_seqs[0] if binder_seqs else ""
            else:
                rec["binder_sequence"] = ":".join(binder_seqs)

            # partner (target/antitarget/self) is mutually exclusive
            partner_list = target_seqs or antitarget_seqs or self_seqs

            if len(partner_list) <= 1:
                rec["target_sequence"] = partner_list[0] if partner_list else ""
            else:
                rec["target_sequence"] = ":".join(partner_list)



            rec.update({
                "binder": binder_dir.name,
                "vs": vs_name,
                "model_idx": int(model_idx),
                "partner": partner_name,
                "target_type": target_type,
            })
            binder_records.append(rec)

    if not binder_records:
        print(f"No valid ipSAE data for {binder_dir.name}")
        return

    df = pd.DataFrame(binder_records)

    csv_path = plots_dir / "ipsae_summary.csv"
    df.to_csv(csv_path, index=False)

    metrics = ["ipSAE_min", "ipSAE_max"]
    partner_order = sorted(df["partner"].dropna().unique().tolist())

    for metric in metrics:
        if metric not in df.columns:
            continue
        plt.figure(figsize=(6, 3.5))
        sns.stripplot(
            data=df,
            x="partner",
            y=metric,
            hue="model_idx",
            alpha=0.7,
            order=partner_order,
        )
        short_title = re.sub(r"^binder_", "", binder_dir.name)
        plt.title(f"{metric} for {short_title}")
        plt.ylabel(metric)
        plt.xlabel("Target / Antitarget / Self")
        plt.xticks(rotation=30)
        handles, labels = plt.gca().get_legend_handles_labels()
        if labels:
            order_idx = sorted(range(len(labels)), key=lambda i: int(labels[i]))
            plt.legend(
                [handles[i] for i in order_idx],
                [labels[i] for i in order_idx],
                title="model_idx",
                loc="best",
            )
            
        plt.tight_layout()
        for ext in ["png"]:
            plt.savefig(plots_dir / f"{metric}_stripplot.{ext}", dpi=200)
        plt.close()

    print(f"Saved: {csv_path}")

import yaml

def extract_sequences_from_yaml(yaml_path):
    with open(yaml_path, "r") as f:
        y = yaml.safe_load(f)

    binder_seqs = []
    target_seqs = []
    antitarget_seqs = []
    self_seqs = []

    for entry in y.get("sequences", []):
        prot = entry.get("protein", {})
        cid = prot.get("id", "")
        seq = prot.get("sequence", "")

        # Binder: A, B, C...
        if len(cid) == 1 and cid.isupper():
            binder_seqs.append(seq)

        # Target: TA, TB, ...
        elif cid.startswith("T") and len(cid) > 1:
            target_seqs.append(seq)

        # Antitarget: AA, AB, ...
        elif cid.startswith("A") and len(cid) > 1:
            antitarget_seqs.append(seq)

        # Self: SA, SB, ...
        elif cid.startswith("S") and len(cid) > 1:
            self_seqs.append(seq)

    return binder_seqs, target_seqs, antitarget_seqs, self_seqs


def plot_overall(root_dir: Path, use_best_model: bool = False):
    """
    Combine all per-binder CSVs and plot heatmaps for ipSAE_min and ipSAE_max.

    All targets & antitargets are pooled together; target_type is kept in the
    DataFrame but not used to split the plots (for now).
    """
    csvs = list(root_dir.glob("binder_*/plots/ipsae_summary.csv"))
    if not csvs:
        print("No binder CSVs found.")
        return

    dfs = []
    for csv in csvs:
        df = pd.read_csv(csv)

        # Backwards compatibility: older CSVs may not have 'partner' or 'target_type'
        if "partner" not in df.columns and "vs" in df.columns:
            df["partner"] = df["vs"].str.extract(r"_vs_(.*)$")
        if "target_type" not in df.columns and "vs" in df.columns:
            df["target_type"] = "unknown"

        df["binder_short"] = df["binder"].str.replace(r"^binder_", "", regex=True)
        dfs.append(df)

    all_df = pd.concat(dfs, ignore_index=True)
    binder_sequences = all_df.groupby("binder_short")["binder_sequence"].first().to_dict()


    # Add the binder_sequence column
    all_df["binder_sequence"] = all_df["binder_short"].map(binder_sequences)

    metrics = [m for m in ["ipSAE_min", "ipSAE_max"] if m in all_df.columns]
    if not metrics:
        print("No ipSAE_min/ipSAE_max metrics found for heatmap plotting.")
        return
    
    # save all_df for reference
    # make dir summary
    Path.mkdir(root_dir/ "summary" ,exist_ok=True)
    all_df_path = root_dir/ "summary" / "ipsae_summary_all_binders.csv"
    all_df.to_csv(all_df_path, index=False)
    print(f"Saved combined ipSAE data at {all_df_path}")

    # ------------------------------------------------------
    # AGGREGATION ACROSS MODELS
    # ------------------------------------------------------
    if use_best_model and "ipSAE_min" in all_df.columns:
        # pick the row (i.e. model) with lowest ipSAE_min per binder/partner
        idx = all_df.groupby(["binder_short", "partner"])["ipSAE_min"].idxmax()
        agg_base = all_df.loc[idx].copy()   # still has target_type
    else:
        # use all rows as base for ordering & averaging
        agg_base = all_df.copy()



    # This is what we actually plot (average across models or best-model rows)
    agg = agg_base.groupby(["binder_short", "partner"])[metrics].mean().reset_index()

    # ------------------------------------------------------
    # ORDER BY BEST BINDING TO TARGET (LOWEST ipSAE_min)
    # ------------------------------------------------------
    # Use only true targets (NOT antitargets) to measure "binding quality"
    targets_only = agg_base[agg_base["target_type"] == "target"]

    # If for some reason no targets exist, fall back to all partners
    source = targets_only if not targets_only.empty else agg_base

    # Binders ordered by highest ipSAE_min (bigger = better)
    binder_order = (
        agg.groupby("binder_short")["ipSAE_min"]
        .max()                        # biggest = best binding
        .sort_values(ascending=False) # best → worst
        .index
        .tolist()
    )

    # Do NOT reorder rows (targets/antitargets)
    type_map = agg_base.groupby("partner")["target_type"].first().to_dict()

    # partner_order = antitargets + targets + self_group
    partner_order = ["self", "target", "antitarget"]

    print("Partner order:", partner_order)

    # Optional: print to verify in logs
    print("Binder order (best→worst by ipSAE_min on targets):", binder_order,"\n")
    # ------------------------------------------------------
    # PLOT HEATMAPS (BOTH USING ipSAE_min-BASED ORDERING)
    # ------------------------------------------------------
    # Define order of partners on Y-axis
    # We want: targets first, then self, then antitargets (or similar logic).
    
    unique_types = ["target", "self", "antitarget"]
    partners_by_type = {t: [] for t in unique_types}
    partners_by_type["unknown"] = []
    
    # Get all partners from agg_base (which has 'partner' and 'target_type')
    unique_partners = agg_base[["partner", "target_type"]].drop_duplicates()

    for _, row in unique_partners.iterrows():
        ptype = row["target_type"]
        pname = row["partner"]
        if ptype in partners_by_type:
            partners_by_type[ptype].append(pname)
        else:
            partners_by_type["unknown"].append(pname)

    # Convert mapping to a flat list in order
    final_partner_order = []
    for t in unique_types + ["unknown"]:
        # sort alphabetically within each group
        final_partner_order.extend(sorted(partners_by_type[t]))

    print(f"Partner order for heatmap: {final_partner_order}")

    for metric in metrics:
        # Use specific 'partner' name instead of class
        # Note: if multiple rows match (binder_short, partner), pivot fails unless aggregated.
        # We already aggregated into 'agg' above by mean().
        
        pivot = agg.pivot(index="partner", columns="binder_short", values=metric)
        
        # Reindex to enforce our sorted order (Target -> Self -> Antitarget)
        # Filter out any partners that might be missing from pivot (shouldn't happen but safe)
        valid_order = [p for p in final_partner_order if p in pivot.index]
        pivot = pivot.reindex(index=valid_order, columns=binder_order)

        plt.figure(figsize=(max(7, len(binder_order) * 0.7),
                            max(5, len(valid_order) * 0.4)))
        sns.heatmap(
            pivot,
            annot=True,
            fmt=".2f",
            cmap="viridis",
            cbar_kws={"label": metric},
            linewidths=0.5,
            vmin=0,
            vmax=1,
        )
        plt.title(metric)
        plt.ylabel("Partner", rotation=90)
        plt.xlabel("Binder")
        plt.yticks(rotation=0)
        # plt.tight_layout()

        for ext in ["png", "svg"]:
            path = root_dir/ "summary"  / f"{metric}_heatmap.{ext}"
            plt.savefig(path, dpi=300, bbox_inches="tight")
            print(f"Saved heatmap for {metric} at {path}")
        plt.close()

        # save CSV (rows = binders, columns = targets)
        csv_out = root_dir/ "summary"  / f"{metric}_heatmap.csv"

        # Add binder_sequence to the final exported heatmap CSV
        # ----------------------------------------------------------
        pivot_out = pivot.T.copy()  # rows = binders

        pivot_out.insert(
            0,
            "binder_sequence",
            pivot_out.index.map(binder_sequences)
        )


        pivot_out.to_csv(csv_out, float_format="%.5f")

        print(f"Saved heatmap data for {metric} as {csv_out}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ipsae_e", type=int, default=15,
                    help="ipSAE PAE cutoff (default: 15 Å)")
    ap.add_argument("--ipsae_d", type=int, default=15,
                    help="ipSAE distance cutoff (default: 15 Å)")
    ap.add_argument("--root_dir", required=True)
    ap.add_argument("--generate_data", action="store_true")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument(
        "--use_best_model",
        action="store_true",
        help="Use only the best model (highest ipSAE_max) per binder/partner instead of averaging"
    )
    ap.add_argument("--num_cpu", type=int, default=1,
                    help="Number of CPUs for parallel processing")
    args = ap.parse_args()


    root = Path(args.root_dir)
    if args.generate_data:
        binder_dirs = [d for d in sorted(root.glob("binder_*")) if d.is_dir()]
        if args.num_cpu == 1:
            # sequential
            for d in binder_dirs:
                analyse_binder(d,args)
        else:
            # parallel
            from multiprocessing import Pool
            with Pool(processes=args.num_cpu) as pool:
                pool.starmap(analyse_binder, [(d, args) for d in binder_dirs])

    if args.plot:
        plot_overall(root, use_best_model=args.use_best_model)


if __name__ == "__main__":
    main()
