"""
Merge all per-method Case 1 prediction CSVs (KEGG-1.8K or Rhea-500) plus the
ground-truth EC labels into one merged_output.csv.

Parametrized (folder/ground-truth/output all via CLI) so it can target
either dataset, and any seed's prediction directory, instead of the
original hardcoded results/Case1/KEGG-1.8K + stale data/KEGG/... path.
"""

import argparse
import glob
import os
from functools import reduce

import pandas as pd


def clean_value(x):
    if pd.isna(x):
        return None
    if isinstance(x, str):
        stripped = x.strip()
        if stripped in ["", "No Significant EC", "No EC Prediction", "No|EC|Prediction", "nan", "nan|", "|nan"]:
            return None
    return x


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge Case 1 per-method prediction CSVs with ground truth.")
    parser.add_argument("--predictions_dir", required=True, help="Folder containing one <Method>.csv per method.")
    parser.add_argument("--ground_truth_csv", required=True, help="CSV/TSV with reaction_id + EC Number columns.")
    parser.add_argument("--methods", nargs="+", default=None,
                         help="Method names to merge, e.g. SIMMER Theia. Defaults to every *.csv in predictions_dir.")
    parser.add_argument("--gt_id_col", default="Reaction ID")
    parser.add_argument("--gt_ec_col", default="EC Number")
    parser.add_argument("--gt_sep", default=",")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.methods:
        files = [os.path.join(args.predictions_dir, f"{m}.csv") for m in args.methods]
        missing = [f for f in files if not os.path.isfile(f)]
        if missing:
            raise SystemExit(f"Method CSV(s) not found: {missing}")
    else:
        files = glob.glob(os.path.join(args.predictions_dir, "*.csv"))
        if not files:
            raise SystemExit(f"No CSV files found in {args.predictions_dir}")

    dataframes = []
    for file in files:
        df = pd.read_csv(file, sep=',', dtype=str)
        first_col = df.columns[0]
        df = df.rename(columns={first_col: "reaction_id"})
        prefix = os.path.splitext(os.path.basename(file))[0] + "_"
        new_columns = {col: prefix + col for col in df.columns if col != "reaction_id"}
        df = df.rename(columns=new_columns)
        dataframes.append(df)

    merged_df = reduce(lambda left, right: pd.merge(left, right, on="reaction_id", how="outer"), dataframes)
    merged_df = merged_df.map(clean_value)

    gt_df = pd.read_csv(args.ground_truth_csv, sep=args.gt_sep, dtype=str)
    gt_df = gt_df.rename(columns={args.gt_id_col: "reaction_id"})
    gt_df = gt_df[["reaction_id", args.gt_ec_col]].rename(columns={args.gt_ec_col: "EC Number"})
    gt_df = gt_df.map(clean_value)

    merged_df = pd.merge(merged_df, gt_df, on="reaction_id", how="outer")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    merged_df.to_csv(args.output, index=False)
    print(f"Merged {len(files)} method CSVs ({len(merged_df)} rows) -> '{args.output}'")


if __name__ == "__main__":
    main()
