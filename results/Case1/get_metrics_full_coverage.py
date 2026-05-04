"""
Evaluate all methods on the subset of reactions where every method produced
at least one prediction (full-coverage evaluation).

This complements get_metrics.py, which penalises missing predictions by
treating them as all-zero binary vectors.
"""

import argparse

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Full-coverage evaluation (reactions predicted by all methods)."
    )
    parser.add_argument(
        "--input", default="results/Case1/merged_output.csv",
        help="Path to merged prediction CSV (default: results/Case1/merged_output.csv)."
    )
    parser.add_argument(
        "--output", default="results/Case1/full_coverage_metrics.csv",
        help="Path for the output evaluation CSV (default: results/Case1/full_coverage_metrics.csv)."
    )
    return parser.parse_args()


def is_invalid_ec_group(s: str) -> bool:
    if not s.strip():
        return True
    for ec in s.split(";")[0].split("|"):
        parts = ec.strip().split(".")
        if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]) and parts[0] != "7":
            return False
    return True


def parse_ecs(s: str) -> set[str]:
    ecs: set[str] = set()
    if not s:
        return ecs
    first_group = s.split(";")[0]
    for ec in first_group.split("|"):
        parts = ec.strip().split(".")
        if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]):
            ecs.add(".".join(parts[:3]))
    return ecs


def evaluate_method(
    col: str,
    df_in: pd.DataFrame,
    mlb: MultiLabelBinarizer,
    all_labels: list[str],
    label_to_class: dict[str, str],
) -> dict:
    y_true = mlb.transform(df_in['true_set'])
    y_pred = mlb.transform(df_in[col + '_pred'])

    support = y_true.sum(axis=0)
    active = support > 0

    mcc_list = []
    for i in range(len(all_labels)):
        if active[i]:
            mcc_list.append(matthews_corrcoef(y_true[:, i], y_pred[:, i]))
        else:
            mcc_list.append(np.nan)

    prec_list = precision_score(y_true, y_pred, average=None, zero_division=0)
    rec_list = recall_score(y_true, y_pred, average=None, zero_division=0)

    active_support = support[active]
    active_mccs = np.array(mcc_list)[active]

    overall_mcc = np.average(active_mccs, weights=active_support)
    overall_prec = np.average(prec_list[active], weights=active_support)
    overall_rec = np.average(rec_list[active], weights=active_support)

    class_mcc: dict[str, float] = {}
    for cls in sorted(set(label_to_class.values())):
        idx = [i for i, lbl in enumerate(all_labels) if label_to_class[lbl] == cls and active[i]]
        if not idx:
            continue
        cls_support = support[idx].sum()
        cls_mccs = [mcc_list[i] for i in idx]
        class_mcc[cls] = float(np.average(cls_mccs, weights=support[idx])) if cls_support > 0 else float('nan')

    return {
        'method': col,
        'n_reactions': len(df_in),
        'overall_mcc': overall_mcc,
        'overall_precision': overall_prec,
        'overall_recall': overall_rec,
        'support': int(support.sum()),
        'class_mcc': class_mcc,
    }


def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.input, dtype=str).fillna("")

    initial_count = len(df)
    invalid_mask = df['EC Number'].apply(is_invalid_ec_group)
    invalid_rows = df[invalid_mask]
    df = df[~invalid_mask].copy()
    filtered_count = len(df)

    print("Eliminated rows with only class 7 or incomplete ECs:")
    print(invalid_rows[['reaction_id', 'EC Number']].to_string(index=False))
    print(f"\nRemoved {initial_count - filtered_count} rows.")
    print(f"Remaining valid reactions: {filtered_count}")

    df['true_set'] = df['EC Number'].apply(parse_ecs)
    prediction_cols = [c for c in df.columns if c not in ['reaction_id', 'EC Number', 'true_set']]

    all_labels: set[str] = set().union(*df['true_set'])
    for col in prediction_cols:
        df[col + '_pred'] = df[col].apply(parse_ecs)
        all_labels |= set().union(*df[col + '_pred'])

    all_labels_sorted = sorted(all_labels)
    mlb = MultiLabelBinarizer(classes=all_labels_sorted)
    mlb.fit(df['true_set'])
    label_to_class = {lbl: lbl.split('.')[0] for lbl in all_labels_sorted}

    print(f"\nReaction counts per method (rows with at least one predicted EC):")
    for col in prediction_cols:
        count = (df[col + '_pred'].apply(lambda x: len(x) > 0)).sum()
        print(f"  {col}: {count} / {filtered_count}")

    all_have_output = pd.Series([True] * len(df), index=df.index)
    for col in prediction_cols:
        all_have_output &= df[col + '_pred'].apply(lambda x: len(x) > 0)

    shared_count = all_have_output.sum()
    print(f"\nReactions where ALL methods produced output: {shared_count} / {filtered_count}")

    df_eval = df[all_have_output].copy()
    print(f"Using these {shared_count} reactions for evaluation.\n")

    results = [evaluate_method(col, df_eval, mlb, all_labels_sorted, label_to_class) for col in prediction_cols]

    for res in results:
        print(f"Method: {res['method']}")
        print(f"  Reactions evaluated:  {res['n_reactions']}")
        print(f"  Weighted MCC:         {res['overall_mcc']:.4f}")
        print(f"  Weighted Precision:   {res['overall_precision']:.4f}")
        print(f"  Weighted Recall:      {res['overall_recall']:.4f}")
        for cls, mcc in res['class_mcc'].items():
            print(f"    Class {cls} MCC:    {mcc:.4f}")
        print()

    records = []
    for res in results:
        for cls, mcc in res['class_mcc'].items():
            records.append({
                'method': res['method'],
                'n_reactions': res['n_reactions'],
                'ec_class': cls,
                'class_mcc': mcc,
                'overall_mcc': res['overall_mcc'],
                'overall_precision': res['overall_precision'],
                'overall_recall': res['overall_recall'],
                'total_support': res['support'],
            })

    pd.DataFrame(records).to_csv(args.output, index=False)
    print(f"Evaluation results saved to '{args.output}'")


if __name__ == "__main__":
    main()
