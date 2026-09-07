"""
Evaluates BEC-Pred's label_assigner.py output against the Rhea Stratified
test split's ground truth: weighted MCC/precision/recall per EC class, class
1 (oxidoreductases) reported first.

BEC-Pred is single-label, so pred_set holds at most one EC per reaction.
"""

import argparse

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer


def parse_pred(ec):
    if not isinstance(ec, str) or not ec.strip():
        return set()
    parts = ec.strip().split(".")
    if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]):
        return {".".join(parts[:3])}
    return set()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", required=True, help="label_assigner.py output CSV (Reaction ID, Prediction)")
    parser.add_argument("--test-tsv", required=True, help="prepared/test.tsv (reaction_id, reaction_smiles, rxn, ec)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    preds = pd.read_csv(args.predictions, dtype=str)
    truth = pd.read_csv(args.test_tsv, sep="\t", dtype=str)

    # Join by ROW POSITION, not by reaction_id value: reaction_id is not unique
    # (multi-label reactions repeat it across rows), and eval_model.py's
    # --reaction_ids_file/--queries contract only guarantees row-order alignment,
    # not id uniqueness -- a value-based merge on reaction_id cartesian-joins
    # duplicate ids (2x2=4 rows instead of 2), silently inflating and
    # double-counting multi-label reactions' support.
    if len(truth) != len(preds):
        raise ValueError(f"truth has {len(truth)} rows but predictions has {len(preds)} -- "
                          "expected row-order-aligned files (test.tsv vs its own queries.txt predictions).")
    assert (truth["reaction_id"].values == preds["reaction_id"].values).all(), \
        "truth/preds reaction_id columns are not row-aligned"
    df = pd.concat([truth.reset_index(drop=True), preds[["Prediction"]].reset_index(drop=True)], axis=1)
    df["true_set"] = df["ec"].apply(lambda ec: {".".join(ec.split(".")[:3])} if isinstance(ec, str) else set())
    df["pred_set"] = df["Prediction"].apply(parse_pred)

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

    class1_mask = df["ec"].astype(str).str.startswith("1.")
    class1_detail_path = args.output.replace(".csv", "_class1_detail.csv")
    df.loc[class1_mask, ["reaction_id", "ec", "Prediction"]].to_csv(class1_detail_path, index=False)
    print(f"Saved class-1 per-reaction detail to {class1_detail_path}")


if __name__ == "__main__":
    main()
