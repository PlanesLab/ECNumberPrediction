"""
Evaluate SelenzymeRF's get_results.py output (Reaction, All ECs) against the
Rhea Stratified test split's ground truth (REACTION_ID, EC_NUMBER -- multiple
rows per REACTION_ID for multi-label reactions), computing the same weighted
MCC/precision/recall metrics used elsewhere in this repo (see
methods/SIMMER/SIMMER_scripts/evaluate_rhea_nocofactor.py), per EC class with
class 1 (oxidoreductases) reported first/separately.
"""

import argparse

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer


def parse_ecs(s):
    """Top-ranked (first ';'-slot) EC(s) from SelenzymeRF's 'All ECs' column."""
    ecs = set()
    if not isinstance(s, str) or not s:
        return ecs
    first_group = s.split(";")[0]
    for ec in first_group.split("|"):
        parts = ec.strip().split(".")
        if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]):
            ecs.add(".".join(parts[:3]))
    return ecs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", required=True, help="get_results.py output CSV (Reaction, All ECs)")
    parser.add_argument("--test-tsv", required=True, help="Rhea Stratified test.tsv (or test_nocofactor.tsv)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    preds = pd.read_csv(args.predictions, dtype=str)
    truth = pd.read_csv(args.test_tsv, sep="\t", dtype=str)

    df = truth.merge(preds, left_on="REACTION_ID", right_on="Reaction", how="left")
    df["true_set"] = df["EC_NUMBER"].apply(lambda ec: {".".join(ec.split(".")[:3])} if isinstance(ec, str) else set())
    df["pred_set"] = df["All ECs"].apply(parse_ecs)

    all_labels = sorted(set().union(*df["true_set"]) | set().union(*df["pred_set"]))
    mlb = MultiLabelBinarizer(classes=all_labels)
    y_true = mlb.fit_transform(df["true_set"])
    y_pred = mlb.transform(df["pred_set"])

    label_to_class = {label: label.split(".")[0] for label in all_labels}
    support = y_true.sum(axis=0)

    mcc_list = [matthews_corrcoef(y_true[:, i], y_pred[:, i]) for i in range(len(all_labels))]
    prec_list = precision_score(y_true, y_pred, average=None, zero_division=0)
    rec_list = recall_score(y_true, y_pred, average=None, zero_division=0)

    overall_mcc = np.average(mcc_list, weights=support)
    overall_prec = np.average(prec_list, weights=support)
    overall_rec = np.average(rec_list, weights=support)
    coverage = (df["pred_set"].apply(len) > 0).mean()

    print(f"Reactions evaluated: {len(df)}")
    print(f"Overall weighted MCC:       {overall_mcc:.4f}")
    print(f"Overall weighted Precision: {overall_prec:.4f}")
    print(f"Overall weighted Recall:    {overall_rec:.4f}")
    print(f"Coverage (>=1 prediction):  {coverage:.2%}")

    records = []
    classes = sorted(set(label_to_class.values()), key=lambda c: (c != "1", c))
    for cls in classes:
        idx = [i for i, lbl in enumerate(all_labels) if label_to_class[lbl] == cls]
        cls_support = int(support[idx].sum())
        cls_mcc = np.average([mcc_list[i] for i in idx], weights=support[idx]) if cls_support > 0 else np.nan
        cls_prec = np.average([prec_list[i] for i in idx], weights=support[idx]) if cls_support > 0 else np.nan
        cls_rec = np.average([rec_list[i] for i in idx], weights=support[idx]) if cls_support > 0 else np.nan
        tag = " <-- FOCUS" if cls == "1" else ""
        print(f"  Class {cls}: MCC={cls_mcc:.4f}  Precision={cls_prec:.4f}  Recall={cls_rec:.4f}  support={cls_support}{tag}")
        records.append({
            "ec_class": cls, "class_mcc": cls_mcc, "class_precision": cls_prec,
            "class_recall": cls_rec, "class_support": cls_support,
            "overall_mcc": overall_mcc, "overall_precision": overall_prec,
            "overall_recall": overall_rec, "coverage": coverage,
            "total_reactions": len(df),
        })

    pd.DataFrame(records).to_csv(args.output, index=False)
    print(f"\nSaved evaluation to {args.output}")

    class1_mask = df["EC_NUMBER"].astype(str).str.startswith("1.")
    class1_detail_path = args.output.replace(".csv", "_class1_detail.csv")
    df.loc[class1_mask, ["REACTION_ID", "EC_NUMBER", "All ECs"]].to_csv(class1_detail_path, index=False)
    print(f"Saved class-1 per-reaction detail to {class1_detail_path}")


if __name__ == "__main__":
    main()
