"""
Script: get_results.py
Author: Josefina Arcagni
Date: 2025-09-11
Description: This script aggregates EC Numbers predicted by SelenzymeRF for each reaction,
             using a combined score from 'Rxn Sim.' and 'Rxn Sim RF.' to rank enzymes.
"""

import os
import glob
import re
import pandas as pd
import argparse

# --- CLI Arguments ---
parser = argparse.ArgumentParser(description="Aggregate EC Numbers by reaction using scoring.")
parser.add_argument("--folder", type=str, required=True, help="Folder containing result CSVs.")
parser.add_argument("--output", type=str, required=True, help="Output CSV file path.")
parser.add_argument("--ec_column", type=str, default="EC Number", help="Column name for EC Numbers.")
args = parser.parse_args()

folder = args.folder
output_csv = args.output
ec_column = args.ec_column

results_all = []

# --- Get all CSV files and sort by KEGG Reaction ID ---
files = glob.glob(os.path.join(folder, "*.csv"))
files = sorted(files, key=lambda x: int(re.search(r"R(\d+)", os.path.basename(x)).group(1)))

for file in files:
    reaction_name = os.path.splitext(os.path.basename(file))[0]  # e.g. R07612
    print(f"Processing: {os.path.basename(file)} → Reaction {reaction_name}")

    df = pd.read_csv(file, header=0)
    df.columns = df.columns.str.strip()  # clean column names

    if ec_column in df.columns:
        # --- Convert numeric columns safely ---
        df["Rxn Sim."] = pd.to_numeric(df.get("Rxn Sim.", 0), errors="coerce").fillna(0)
        df["Rxn Sim RF."] = pd.to_numeric(df.get("Rxn Sim RF.", 0), errors="coerce").fillna(0)

        # --- Combined scoring rule ---
        def combined_score(row):
            sim, sim_rf = row["Rxn Sim."], row["Rxn Sim RF."]
            if sim > 0 and sim_rf > 0:
                return 0.5 * sim + 0.5 * sim_rf
            elif sim > 0:
                return sim
            elif sim_rf > 0:
                return sim_rf
            else:
                return 0.0

        df["Combined_Score"] = df.apply(combined_score, axis=1)

        # --- Sort by score ---
        df_sorted = df.sort_values(by="Combined_Score", ascending=False)

        # --- Clean EC numbers ---
        def process_ec(ec):
            if pd.isna(ec) or str(ec).strip() in ["", "-"]:
                return ""
            ecs = [e.strip() for e in str(ec).split(";") if e.strip()]
            return "|".join(sorted(set(ecs)))

        df_sorted[ec_column] = df_sorted[ec_column].apply(process_ec)

        # --- Collect unique ECs (ranked order preserved) ---
        unique_ecs = [ec for ec in df_sorted[ec_column].unique() if ec]
        all_ecs = ";".join(unique_ecs)
    else:
        print(f"⚠️ No EC Number column in {file}")
        all_ecs = ""

    results_all.append({"Reaction": reaction_name, "All ECs": all_ecs})

# --- Save final aggregated results ---
df_all = pd.DataFrame(results_all)
df_all.to_csv(output_csv, index=False)

print("\n✅ All ECs by Reaction saved to:", output_csv)
print(df_all.head(10))
