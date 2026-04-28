"""
Compute weighted multilabel MCC for BEC-Pred, Theia (top1), and Gemma4EC
across the three Rhea splits (Scaffold, Stratified, Time).
Uses the same methodology as get_metrics.py.
"""
import pandas as pd
import numpy as np
from sklearn.metrics import matthews_corrcoef, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer
import warnings
warnings.filterwarnings("ignore")

SPLITS = {
    "Scaffold": {
        "test_tsv":   "/scratch/jarcagniriv/Case1-Rhea/Database/Scaffold/test.tsv",
        "becpred":    "/scratch/jarcagniriv/ECNumberPrediction/results/Case2/results-splits/Scaffold_becpred.csv",
        "theia":      "/scratch/jarcagniriv/Case1/theia/results/scaffold_ec123_test_predictions.csv",
        "gemma4ec":   "/scratch/jarcagniriv/gemma4ec/Rhea_Scaffold/results/predictions_pred.EC_NUMBER.EC_NUMBER_PROMPT_2.20260416-212127.tsv",
    },
    "Stratified": {
        "test_tsv":   "/scratch/jarcagniriv/Case1-Rhea/Database/Stratified/test.tsv",
        "becpred":    "/scratch/jarcagniriv/ECNumberPrediction/results/Case2/results-splits/Stratified_becpred.csv",
        "theia":      "/scratch/jarcagniriv/Case1/theia/results/stratified_ec123_test_predictions.csv",
        "gemma4ec":   "/scratch/jarcagniriv/gemma4ec/Rhea_Stratified/results/predictions_pred.EC_NUMBER.EC_NUMBER_PROMPT_2.20260416-213600.tsv",
    },
    "Time": {
        "test_tsv":   "/scratch/jarcagniriv/Case1-Rhea/Database/Time/test.tsv",
        "becpred":    "/scratch/jarcagniriv/ECNumberPrediction/results/Case2/results-splits/Time_becpred.csv",
        "theia":      "/scratch/jarcagniriv/Case1/theia/results/time_ec123_test_predictions.csv",
        "gemma4ec":   "/scratch/jarcagniriv/gemma4ec/Rhea_Time/results/predictions_pred.EC_NUMBER.EC_NUMBER_PROMPT_2.20260416-214016.tsv",
    },
}


def parse_ecs(s):
    ecs = set()
    if not s or (isinstance(s, float) and np.isnan(s)):
        return ecs
    for group in str(s).split(";"):
        for ec in group.split("|"):
            parts = ec.strip().split(".")
            if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]):
                norm = ".".join(str(int(p)) for p in parts[:3])
                if norm.split(".")[0] != "0":
                    ecs.add(norm)
    return ecs


def evaluate_methods(df, prediction_cols):
    all_labels = set().union(*df["true_set"])
    for col in prediction_cols:
        all_labels |= set().union(*df[col + "_pred"])
    all_labels = sorted(all_labels)

    if not all_labels:
        return {}

    mlb = MultiLabelBinarizer(classes=all_labels)
    mlb.fit(df["true_set"])
    label_to_class = {lbl: lbl.split(".")[0] for lbl in all_labels}

    results = {}
    for col in prediction_cols:
        y_true = mlb.transform(df["true_set"])
        y_pred = mlb.transform(df[col + "_pred"])
        support = y_true.sum(axis=0)

        mcc_list = [matthews_corrcoef(y_true[:, i], y_pred[:, i]) for i in range(len(all_labels))]
        prec_list = precision_score(y_true, y_pred, average=None, zero_division=0)
        rec_list = recall_score(y_true, y_pred, average=None, zero_division=0)

        overall_mcc = float(np.average(mcc_list, weights=support))
        overall_prec = float(np.average(prec_list, weights=support))
        overall_rec = float(np.average(rec_list, weights=support))
        coverage = float((df[col + "_pred"].apply(len) > 0).mean())

        class_mcc = {}
        for cls in sorted(set(label_to_class.values())):
            idx = [i for i, lbl in enumerate(all_labels) if label_to_class[lbl] == cls]
            if not idx:
                continue
            cls_support = support[idx].sum()
            cls_mccs = [mcc_list[i] for i in idx]
            class_mcc[cls] = float(np.average(cls_mccs, weights=support[idx])) if cls_support > 0 else float("nan")

        results[col] = {
            "overall_mcc": overall_mcc,
            "overall_precision": overall_prec,
            "overall_recall": overall_rec,
            "coverage": coverage,
            "total_support": int(support.sum()),
            "class_mcc": class_mcc,
        }
    return results


all_records = []

for split, paths in SPLITS.items():
    print(f"\n=== {split} ===")

    # Ground truth
    test = pd.read_csv(paths["test_tsv"], sep="\t")
    test.columns = [c.strip() for c in test.columns]
    test = test.rename(columns={"REACTION_ID": "reaction_id", "REACTION_SMILES": "smiles", "EC_NUMBER": "ec_number"})
    test["true_set"] = test["ec_number"].apply(parse_ecs)
    test = test[test["true_set"].apply(len) > 0].copy()

    # BEC-Pred (join by reaction_id)
    bec = pd.read_csv(paths["becpred"])
    bec.columns = [c.strip() for c in bec.columns]
    bec = bec.rename(columns={"Reaction ID": "reaction_id", "Prediction": "BEC-Pred"})
    bec["reaction_id"] = bec["reaction_id"].astype(str)
    test["reaction_id"] = test["reaction_id"].astype(str)
    df = test.merge(bec[["reaction_id", "BEC-Pred"]], on="reaction_id", how="left")

    # Theia (join by smiles)
    theia = pd.read_csv(paths["theia"])
    theia = theia[["rxn_smiles", "top1_ec"]].rename(columns={"rxn_smiles": "smiles", "top1_ec": "Theia"})
    df = df.merge(theia, on="smiles", how="left")

    # Gemma4EC (join by smiles)
    gemma = pd.read_csv(paths["gemma4ec"], sep="\t")
    gemma = gemma[["smiles", "pred_fragment"]].rename(columns={"pred_fragment": "Gemma4EC"})
    df = df.merge(gemma, on="smiles", how="left")

    # Fill missing with empty string (penalizes no prediction)
    for col in ["BEC-Pred", "Theia", "Gemma4EC"]:
        df[col] = df[col].fillna("")
        df[col + "_pred"] = df[col].apply(parse_ecs)

    n = len(df)
    for col in ["BEC-Pred", "Theia", "Gemma4EC"]:
        covered = (df[col + "_pred"].apply(len) > 0).sum()
        print(f"  {col}: {covered}/{n} reactions with predictions")

    results = evaluate_methods(df, ["BEC-Pred", "Theia", "Gemma4EC"])

    for method, res in results.items():
        print(f"  {method}: MCC={res['overall_mcc']:.4f}, Prec={res['overall_precision']:.4f}, Rec={res['overall_recall']:.4f}, Cov={res['coverage']:.4f}")
        for cls, mcc in res["class_mcc"].items():
            all_records.append({
                "split": split,
                "method": method,
                "ec_class": cls,
                "class_mcc": mcc,
                "overall_mcc": res["overall_mcc"],
                "overall_precision": res["overall_precision"],
                "overall_recall": res["overall_recall"],
                "coverage": res["coverage"],
                "total_support": res["total_support"],
            })

out = pd.DataFrame(all_records)
out_path = "/scratch/jarcagniriv/ECNumberPrediction/results/Case2/evaluation_summary_splits.csv"
out.to_csv(out_path, index=False)
print(f"\nSaved to {out_path}")

print("\n=== OVERALL MCC SUMMARY ===")
summary = out.drop_duplicates(subset=["split", "method"])[["split", "method", "overall_mcc", "coverage", "total_support"]]
print(summary.to_string(index=False))
