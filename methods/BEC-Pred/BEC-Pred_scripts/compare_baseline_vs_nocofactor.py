"""
Build a baseline-vs-no-cofactor comparison CSV from two
evaluate_rhea_stratified.py outputs, in the same format as
results/Case2/SIMMER_baseline_vs_nocofactor_comparison.csv.
"""

import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--nocofactor", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    base = pd.read_csv(args.baseline)
    nocof = pd.read_csv(args.nocofactor)

    merged = base[["ec_class", "class_mcc", "class_precision", "class_recall", "class_support"]].merge(
        nocof[["ec_class", "class_mcc", "class_precision", "class_recall", "class_support"]],
        on="ec_class", suffixes=("_baseline", "_nocofactor"),
    )
    merged = merged.rename(columns={
        "class_mcc_baseline": "baseline_mcc", "class_mcc_nocofactor": "nocofactor_mcc",
        "class_precision_baseline": "baseline_precision", "class_precision_nocofactor": "nocofactor_precision",
        "class_recall_baseline": "baseline_recall", "class_recall_nocofactor": "nocofactor_recall",
        "class_support_baseline": "baseline_support", "class_support_nocofactor": "nocofactor_support",
    })
    merged["mcc_delta"] = (merged["nocofactor_mcc"] - merged["baseline_mcc"]).round(4)
    merged["precision_delta"] = (merged["nocofactor_precision"] - merged["baseline_precision"]).round(4)
    merged["recall_delta"] = (merged["nocofactor_recall"] - merged["baseline_recall"]).round(4)

    cols = ["ec_class", "baseline_mcc", "nocofactor_mcc", "mcc_delta",
            "baseline_precision", "nocofactor_precision", "precision_delta",
            "baseline_recall", "nocofactor_recall", "recall_delta",
            "baseline_support", "nocofactor_support"]
    merged = merged[cols]

    overall = {
        "ec_class": "OVERALL",
        "baseline_mcc": base["overall_mcc"].iloc[0], "nocofactor_mcc": nocof["overall_mcc"].iloc[0],
        "baseline_precision": base["overall_precision"].iloc[0], "nocofactor_precision": nocof["overall_precision"].iloc[0],
        "baseline_recall": base["overall_recall"].iloc[0], "nocofactor_recall": nocof["overall_recall"].iloc[0],
        "baseline_support": base["total_reactions"].iloc[0], "nocofactor_support": nocof["total_reactions"].iloc[0],
    }
    overall["mcc_delta"] = round(overall["nocofactor_mcc"] - overall["baseline_mcc"], 4)
    overall["precision_delta"] = round(overall["nocofactor_precision"] - overall["baseline_precision"], 4)
    overall["recall_delta"] = round(overall["nocofactor_recall"] - overall["baseline_recall"], 4)

    out = pd.concat([merged, pd.DataFrame([overall])[cols]], ignore_index=True)
    out.to_csv(args.output, index=False)
    print(out.to_string(index=False))
    print(f"\nSaved to {args.output}")


if __name__ == "__main__":
    main()
