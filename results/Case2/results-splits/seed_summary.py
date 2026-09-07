"""
Aggregate results/Case2/results-splits/seed_metrics_{Scaffold,Stratified}.csv (each written by
score_seed.py, one row per method x seed) into a single mean/std-over-seeds summary.

Companion to results/Case2/bootstrap_summary.csv's aggregation logic, but over genuinely
different reruns (3 real train/test splits) rather than statistical resampling of one fixed
merged_output.csv -- same idea as results/Case1's seed-run summary.
"""

import pandas as pd

SPLITS = ["Scaffold", "Stratified"]
INPUT_TEMPLATE = "results/Case2/results-splits/seed_metrics_{split}.csv"
OUTPUT = "results/Case2/results-splits/seed_summary.csv"


def main() -> None:
    frames = []
    for split in SPLITS:
        path = INPUT_TEMPLATE.format(split=split)
        try:
            df = pd.read_csv(path)
        except FileNotFoundError:
            print(f"Skipping '{split}': '{path}' not found yet (run run_join_and_score.sh first)")
            continue
        df["split"] = split
        frames.append(df)

    if not frames:
        raise SystemExit("No seed_metrics_*.csv files found -- nothing to aggregate.")

    all_df = pd.concat(frames, ignore_index=True)
    n_seeds = all_df.groupby(["split", "method"])["seed"].nunique()
    incomplete = n_seeds[n_seeds < 3]
    if len(incomplete) > 0:
        print("NOTE: these method/split combos have fewer than 3 seeds so far:")
        print(incomplete.to_string())

    summary = all_df.groupby(["split", "method"])[["mcc", "ppv", "recall"]].agg(["mean", "std"])
    summary.columns = ["_".join(c) for c in summary.columns]
    summary = summary.reset_index()
    summary.to_csv(OUTPUT, index=False)
    print(f"Wrote '{OUTPUT}'")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
