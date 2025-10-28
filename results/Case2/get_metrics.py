import pandas as pd
import numpy as np
from sklearn.metrics import matthews_corrcoef, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer
from datetime import datetime

# === File path ===
csv_path = 'results/Case2/merged_output.csv'

# === Load & Clean Data ===
df = pd.read_csv(csv_path, dtype=str).fillna("")

def is_invalid_ec_group(s):
    if not s.strip():
        return True
    for ec in s.split(";")[0].split("|"):
        parts = ec.strip().split(".")
        if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]) and parts[0] != "7":
            return False
    return True

initial_count = len(df)
invalid_mask = df['EC Number'].apply(is_invalid_ec_group)
invalid_rows = df[invalid_mask]
df = df[~invalid_mask].copy()
filtered_count = len(df)

print("❌ Eliminated rows with only class 7 or incomplete ECs:")
print(invalid_rows[['reaction_id', 'EC Number']].to_string(index=False))
print(f"\n🔎 Removed {initial_count - filtered_count} rows.")
print(f"✅ Remaining valid reactions: {filtered_count}")

# === Parse ECs ===
def parse_ecs(s):
    ecs = set()
    if not s:
        return ecs
    first_group = s.split(";")[0]
    for ec in first_group.split("|"):
        parts = ec.strip().split(".")
        if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]):
            ecs.add(".".join(parts[:3]))
    return ecs

df['true_set'] = df['EC Number'].apply(parse_ecs)
prediction_cols = [c for c in df.columns if c not in ['reaction_id', 'EC Number', 'true_set']]

# === Collect All Labels ===
all_labels = set().union(*df['true_set'])
for col in prediction_cols:
    df[col + '_pred'] = df[col].apply(parse_ecs)
    all_labels |= set().union(*df[col + '_pred'])

all_labels = sorted(all_labels)
mlb = MultiLabelBinarizer(classes=all_labels)
mlb.fit(df['true_set'])

label_to_class = {label: label.split('.')[0] for label in all_labels}

# === Evaluation Function ===
def evaluate_method(col):
    y_true = mlb.transform(df['true_set'])
    
    # Penalize missing predictions (empty set means all zeros)
    y_pred_raw = df[col + '_pred'].apply(lambda x: x if x else set())
    y_pred = mlb.transform(y_pred_raw)

    mcc_list = [matthews_corrcoef(y_true[:, i], y_pred[:, i]) for i in range(len(all_labels))]
    prec_list = precision_score(y_true, y_pred, average=None, zero_division=0)
    rec_list = recall_score(y_true, y_pred, average=None, zero_division=0)
    support = y_true.sum(axis=0)

    overall_mcc = np.average(mcc_list, weights=support)
    overall_prec = np.average(prec_list, weights=support)
    overall_rec = np.average(rec_list, weights=support)

    class_mcc = {}
    for cls in sorted(set(label_to_class.values())):
        idx = [i for i, lbl in enumerate(all_labels) if label_to_class[lbl] == cls]
        if not idx:
            continue
        cls_support = support[idx].sum()
        cls_mccs = [mcc_list[i] for i in idx]
        class_mcc[cls] = np.average(cls_mccs, weights=support[idx]) if cls_support > 0 else np.nan

    # === Coverage: rows with at least one prediction ===
    coverage = (df[col + '_pred'].apply(lambda x: len(x) > 0)).mean()

    return {
        'method': col,
        'overall_mcc': overall_mcc,
        'overall_precision': overall_prec,
        'overall_recall': overall_rec,
        'coverage': coverage,
        'support': int(support.sum()),
        'class_mcc': class_mcc
    }

# === Run Evaluation ===
results = [evaluate_method(col) for col in prediction_cols]

# === Display Results ===
for res in results:
    print(f"\nMethod: {res['method']}")
    print(f"  • Weighted MCC:       {res['overall_mcc']:.4f}")
    print(f"  • Weighted Precision: {res['overall_precision']:.4f}")
    print(f"  • Weighted Recall:    {res['overall_recall']:.4f}")
    print(f"  • Coverage:           {res['coverage']:.2%}")
    for cls, mcc in res['class_mcc'].items():
        print(f"    – Class {cls} MCC:  {mcc:.4f}")

# === Export to CSV ===
records = []
for res in results:
    for cls, mcc in res['class_mcc'].items():
        records.append({
            'method': res['method'],
            'ec_class': cls,
            'class_mcc': mcc,
            'overall_mcc': res['overall_mcc'],
            'overall_precision': res['overall_precision'],
            'overall_recall': res['overall_recall'],
            'coverage': res['coverage'],
            'total_support': res['support']
        })

results_df = pd.DataFrame(records)
filename = "results/Case2/evaluation_summary.csv"
results_df.to_csv(filename, index=False)
print(f"\n📁 Evaluation results saved to '{filename}'")