"""
Merge all per-method Case 2 (MetaNetX) prediction CSVs plus the ground-truth
EC labels into one merged_output.csv.

Parametrized (folder/ground-truth/output all via CLI) so it can target any
seed's prediction directory, instead of the original hardcoded
results/Case2/results + stale data/MetaNetX/test_reactions.tsv path.
Same merge logic as results/Case1/join_results.py, just defaulted for
Case 2's tab-separated ground truth (reaction_id/ec columns).
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
    parser = argparse.ArgumentParser(description="Merge Case 2 per-method prediction CSVs with ground truth.")
    parser.add_argument("--predictions_dir", required=True, help="Folder containing one <Method>.csv per method.")
    parser.add_argument("--ground_truth_csv", required=True, help="TSV with reaction_id + ec columns (e.g. a seed's test.tsv).")
    parser.add_argument("--methods", nargs="+", default=None,
                         help="Method names to merge, e.g. SIMMER Theia. Defaults to every *.csv in predictions_dir.")
    parser.add_argument("--gt_id_col", default="reaction_id")
    parser.add_argument("--gt_ec_col", default="ec")
    parser.add_argument("--gt_sep", default="\t")
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
        # Reactions with multiple EC labels get one row per EC in the source data, so the same
        # reaction_id can appear more than once in a prediction file. A method's prediction is
        # driven only by the reaction SMILES (identical across those duplicate rows), so keeping
        # more than one is redundant -- worse, if two DIFFERENT prediction files both have
        # duplicate rows for the same reaction_id, pd.merge's outer join builds the full
        # cross-product of them (17 x 17 rows for a 17-EC reaction), which is both wrong (massively
        # over-weights that reaction) and can blow up runtime/memory. Dedup each file to one row
        # per reaction_id before merging; ground truth's duplicate EC rows are untouched, so each
        # still gets scored against the (now singular, correct) prediction.
        df = df.drop_duplicates(subset="reaction_id", keep="first")
        prefix = os.path.splitext(os.path.basename(file))[0] + "_"
        new_columns = {col: prefix + col for col in df.columns if col != "reaction_id"}
        df = df.rename(columns=new_columns)
        dataframes.append(df)

    merged_df = reduce(lambda left, right: pd.merge(left, right, on="reaction_id", how="outer"), dataframes)
    merged_df = merged_df.map(clean_value)

    gt_df = pd.read_csv(args.ground_truth_csv, sep=args.gt_sep, dtype=str)
    gt_df = gt_df.rename(columns={args.gt_id_col: "reaction_id"})
    gt_df = gt_df[["reaction_id", args.gt_ec_col]].rename(columns={args.gt_ec_col: "ec"})
    gt_df = gt_df.map(clean_value)

    merged_df = pd.merge(merged_df, gt_df, on="reaction_id", how="outer")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    merged_df.to_csv(args.output, index=False)
    print(f"Merged {len(files)} method CSVs ({len(merged_df)} rows) -> '{args.output}'")


if __name__ == "__main__":
    main()
