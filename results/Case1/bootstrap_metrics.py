"""
Bootstrap-resampled MCC/PPV/Recall for Case 1 (KEGG-1.8K).

For each of --n-seeds seeds, resamples the (already-filtered) reaction rows
with replacement and recomputes weighted MCC/precision/recall per method.
Label space (all_labels) is fixed from the full filtered dataset so every
resample is scored against the same class set.
"""

import argparse

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer

PREDICTION_COLS = [
    'E-zyme1', 'E-zyme2', 'BridgIT', 'SelenzymeRF', 'SIMMER', 'Theia', 'BEC-Pred',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Bootstrap MCC/PPV/Recall for Case 1.")
    parser.add_argument("--input", default="results/Case1/merged_output.csv")
    parser.add_argument("--output", default="results/Case1/bootstrap_metrics.csv")
    parser.add_argument("--summary-output", default="results/Case1/bootstrap_summary.csv")
    parser.add_argument("--n-seeds", type=int, default=10)
    return parser.parse_args()


def is_invalid_ec_group(s: str) -> bool:
    if not s.strip():
        return True
    for ec in s.split(";")[0].split("|"):
        parts = ec.strip().split(".")
        if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]) and parts[0] != "7":
            return False
    return True


def parse_ecs(s: str) -> set:
    ecs = set()
    if not s:
        return ecs
    for ec in s.split(";")[0].split("|"):
        parts = ec.strip().split(".")
        if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]):
            ecs.add(".".join(parts[:3]))
    return ecs


def weighted_metrics(y_true, y_pred, n_labels) -> tuple:
    mcc_list = [matthews_corrcoef(y_true[:, i], y_pred[:, i]) for i in range(n_labels)]
    prec_list = precision_score(y_true, y_pred, average=None, zero_division=0)
    rec_list = recall_score(y_true, y_pred, average=None, zero_division=0)
    support = y_true.sum(axis=0)
    if support.sum() == 0:
        return float('nan'), float('nan'), float('nan')
    mcc = np.average(mcc_list, weights=support)
    ppv = np.average(prec_list, weights=support)
    rec = np.average(rec_list, weights=support)
    return mcc, ppv, rec


def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.input, dtype=str).fillna("")
    df = df[~df['EC Number'].apply(is_invalid_ec_group)].copy().reset_index(drop=True)

    df['true_set'] = df['EC Number'].apply(parse_ecs)
    all_labels = set().union(*df['true_set'])
    for col in PREDICTION_COLS:
        df[col + '_pred'] = df[col].apply(parse_ecs)
        all_labels |= set().union(*df[col + '_pred'])
    all_labels = sorted(all_labels)

    mlb = MultiLabelBinarizer(classes=all_labels)
    mlb.fit(df['true_set'])
    n_labels = len(all_labels)
    n = len(df)

    y_true_full = mlb.transform(df['true_set'])
    y_pred_full = {col: mlb.transform(df[col + '_pred']) for col in PREDICTION_COLS}

    print(f"Filtered reactions: {n}, labels: {n_labels}, seeds: {args.n_seeds}")

    records = []
    for seed in range(args.n_seeds):
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, n, size=n)  # resample WITH replacement, shared across methods
        y_true = y_true_full[idx]
        for col in PREDICTION_COLS:
            y_pred = y_pred_full[col][idx]
            mcc, ppv, rec = weighted_metrics(y_true, y_pred, n_labels)
            records.append({'method': col, 'seed': seed, 'mcc': mcc, 'ppv': ppv, 'recall': rec})
            print(f"  seed={seed} {col:12s} MCC={mcc:.4f} PPV={ppv:.4f} Recall={rec:.4f}")

    results_df = pd.DataFrame(records)
    results_df.to_csv(args.output, index=False)
    print(f"\nPer-seed results saved to '{args.output}'")

    summary = results_df.groupby('method')[['mcc', 'ppv', 'recall']].agg(['mean', 'std'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    summary = summary.reset_index()
    summary.to_csv(args.summary_output, index=False)
    print(f"Summary (mean/std over {args.n_seeds} seeds) saved to '{args.summary_output}'")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
