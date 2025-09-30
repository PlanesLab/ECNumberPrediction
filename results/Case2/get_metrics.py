import pandas as pd
import numpy as np
from sklearn.metrics import matthews_corrcoef, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer
import math
import warnings

# --- Configuration ---
CSV_PATH = '/scratch/jarcagniriv/Case2/Results/results/merged_output.csv'
OUT_SUMMARY = '/scratch/jarcagniriv/Case2/Results/evaluation_summary.csv'
OUT_MISSING = '/scratch/jarcagniriv/Case2/Results/missing_coverage_per_method.csv'

# --- Utilities for EC parsing & cleaning ---
def is_invalid_ec_group(s: str) -> bool:
    if not s or not s.strip():
        return True
    first_group = s.split(";")[0]
    for ec in first_group.split("|"):
        parts = [p.strip() for p in ec.strip().split(".")]
        if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]) and parts[0] != "7":
            return False
    return True

def parse_ecs_three_part(s: str) -> set:
    ecs = set()
    if not s or not s.strip():
        return ecs
    first_group = s.split(";")[0]
    for ec in first_group.split("|"):
        parts = [p.strip() for p in ec.strip().split(".")]
        if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]):
            ecs.add(".".join(parts[:3]))
    return ecs

# --- Top-1 helpers (levels) ---
def get_top1_group_str(s: str) -> str:
    if not s or not s.strip():
        return ""
    return s.split(";")[0].strip()

def split_group_to_ecs(group_str: str) -> list:
    if not group_str:
        return []
    return [ec.strip() for ec in group_str.split("|") if ec.strip()]

def truncate_ec_to_level(ec: str, level: int):
    if not ec:
        return None
    parts = [p.strip() for p in ec.split(".")]
    if len(parts) < level:
        return None
    if not all(p.isdigit() for p in parts[:level]):
        return None
    return ".".join(parts[:level])

# --- Top-1 MCC computation ---
def compute_top1_mcc_weighted_by_label(df_local: pd.DataFrame, pred_col: str):
    """
    Compute weighted Top-1 MCC for levels 1..4.
    Each truncated EC (level 1..4) is treated as a binary label, MCC computed per-label,
    and then averaged weighted by label support.
    """
    results = {}
    for level in [1, 2, 3, 4]:
        y_true_sets = []
        y_pred_sets = []
        for _, row in df_local.iterrows():
            true_group = get_top1_group_str(row['ec'])
            true_ecs = split_group_to_ecs(true_group)
            true_trunc = {truncate_ec_to_level(ec, level) for ec in true_ecs if truncate_ec_to_level(ec, level)}
            y_true_sets.append(true_trunc)

            pred_group = get_top1_group_str(row.get(pred_col, ""))
            pred_ecs = split_group_to_ecs(pred_group)
            pred_trunc = {truncate_ec_to_level(ec, level) for ec in pred_ecs if truncate_ec_to_level(ec, level)}
            y_pred_sets.append(pred_trunc)

        all_level_labels = set().union(*y_true_sets, *y_pred_sets) if len(y_true_sets) > 0 else set()
        if not all_level_labels:
            results[f"mcc_level{level}"] = np.nan
            continue

        all_level_labels = sorted(all_level_labels)
        mlb_lvl = MultiLabelBinarizer(classes=all_level_labels)
        mlb_lvl.fit(y_true_sets + y_pred_sets)
        Y_true = mlb_lvl.transform(y_true_sets)
        Y_pred = mlb_lvl.transform(y_pred_sets)

        mcc_list = []
        support = Y_true.sum(axis=0)

        for i in range(len(all_level_labels)):
            y_t = Y_true[:, i]
            y_p = Y_pred[:, i]
            try:
                if (y_t.sum() == 0 and y_p.sum() == 0) or (y_t.sum() == len(y_t) and y_p.sum() == len(y_p)):
                    mcc_val = 0.0
                else:
                    mcc_val = matthews_corrcoef(y_t, y_p)
                    if isinstance(mcc_val, float) and math.isnan(mcc_val):
                        mcc_val = 0.0
            except Exception:
                mcc_val = 0.0
            mcc_list.append(mcc_val)

        if support.sum() > 0:
            weighted_mcc = float(np.average(mcc_list, weights=support))
        else:
            weighted_mcc = np.nan

        results[f"mcc_level{level}"] = weighted_mcc

    return results

# --- Main evaluation (unchanged general MCC) ---
def evaluate_method(df_eval: pd.DataFrame, col: str, all_labels: list, label_to_class: dict):
    if all_labels:
        mlb = MultiLabelBinarizer(classes=all_labels)
        mlb.fit([])
        Y_true = mlb.transform(df_eval['true_set'])
    else:
        Y_true = np.zeros((len(df_eval), 0), dtype=int)

    y_pred_raw = df_eval[col + '_pred'].apply(lambda x: x if x else set())
    if all_labels:
        Y_pred = mlb.transform(y_pred_raw)
    else:
        Y_pred = np.zeros((len(df_eval), 0), dtype=int)

    n_labels = Y_true.shape[1]
    mcc_list = []
    support = np.array([])

    if n_labels > 0:
        for i in range(n_labels):
            y_t = Y_true[:, i]
            y_p = Y_pred[:, i]
            try:
                if (y_t.sum() == 0 and y_p.sum() == 0) or (y_t.sum() == len(y_t) and y_p.sum() == len(y_p)):
                    mcc_val = 0.0
                else:
                    mcc_val = matthews_corrcoef(y_t, y_p)
                    if isinstance(mcc_val, float) and math.isnan(mcc_val):
                        mcc_val = 0.0
            except Exception:
                mcc_val = 0.0
            mcc_list.append(mcc_val)

        prec_list = precision_score(Y_true, Y_pred, average=None, zero_division=0)
        rec_list = recall_score(Y_true, Y_pred, average=None, zero_division=0)
        support = Y_true.sum(axis=0)
    else:
        prec_list = np.array([])
        rec_list = np.array([])

    total_support = int(support.sum()) if support.size else 0
    if support.size and support.sum() > 0:
        overall_mcc = float(np.average(mcc_list, weights=support))
        overall_prec = float(np.average(prec_list, weights=support))
        overall_rec = float(np.average(rec_list, weights=support))
    else:
        overall_mcc = np.nan
        overall_prec = np.nan
        overall_rec = np.nan

    class_mcc = {}
    for cls in sorted(set(label_to_class.values())):
        idx = [i for i, lbl in enumerate(all_labels) if label_to_class[lbl] == cls]
        if not idx:
            continue
        cls_support = support[idx].sum() if support.size else 0
        if cls_support > 0:
            cls_mccs = [mcc_list[i] for i in idx]
            class_mcc[cls] = float(np.average(cls_mccs, weights=support[idx]))
        else:
            class_mcc[cls] = np.nan

    coverage = (df_eval[col + '_pred'].apply(lambda x: len(x) > 0)).mean()

    # Top-1 MCC (new)
    top1_mcc = compute_top1_mcc_weighted_by_label(df_eval, col)

    return {
        'method': col,
        'overall_mcc': overall_mcc,
        'overall_precision': overall_prec,
        'overall_recall': overall_rec,
        'coverage': coverage,
        'support': total_support,
        'class_mcc': class_mcc,
        **top1_mcc
    }
