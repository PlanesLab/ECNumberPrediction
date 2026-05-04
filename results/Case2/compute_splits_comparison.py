"""
Compute weighted multilabel MCC for BEC-Pred, Theia, and Gemma4EC across
the three Rhea splits (Scaffold, Stratified, Time).

Uses the same methodology as get_metrics.py.

External result files (Theia, Gemma4EC) live outside this repository and
must be provided via CLI arguments. BEC-Pred results are read from this
repo's results/Case2/results-splits/ directory.

Example usage:
    python compute_splits_comparison.py \
        --rhea_data_dir  /path/to/rhea/splits \
        --theia_dir      /path/to/theia/results \
        --gemma_dir      /path/to/gemma4ec/results \
        --output         results/Case2/evaluation_summary_splits.csv
"""

import argparse
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer

warnings.filterwarnings("ignore")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare BEC-Pred, Theia, and Gemma4EC across Rhea splits."
    )
    parser.add_argument(
        "--rhea_data_dir", required=True,
        help=(
            "Base directory containing the Rhea split sub-folders "
            "(Scaffold/, Stratified/, Time/), each with a test.tsv."
        ),
    )
    parser.add_argument(
        "--theia_dir", required=True,
        help=(
            "Directory containing Theia prediction CSVs named "
            "{scaffold,stratified,time}_ec123_test_predictions.csv."
        ),
    )
    parser.add_argument(
        "--gemma_dir", required=True,
        help=(
            "Base directory for Gemma4EC results. Expects sub-folders "
            "Rhea_Scaffold/, Rhea_Stratified/, Rhea_Time/ each containing "
            "a predictions TSV matching 'predictions_pred.EC_NUMBER.*.tsv'."
        ),
    )
    parser.add_argument(
        "--becpred_dir",
        default="results/Case2/results-splits",
        help=(
            "Directory containing BEC-Pred CSVs named "
            "{Scaffold,Stratified,Time}_becpred.csv "
            "(default: results/Case2/results-splits)."
        ),
    )
    parser.add_argument(
        "--output",
        default="results/Case2/evaluation_summary_splits.csv",
        help="Path for the output CSV (default: results/Case2/evaluation_summary_splits.csv).",
    )
    return parser.parse_args()


def parse_ecs(s) -> set[str]:
    ecs: set[str] = set()
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


def evaluate_methods(df: pd.DataFrame, prediction_cols: list[str]) -> dict:
    all_labels: set[str] = set().union(*df["true_set"])
    for col in prediction_cols:
        all_labels |= set().union(*df[col + "_pred"])
    all_labels_sorted = sorted(all_labels)

    if not all_labels_sorted:
        return {}

    mlb = MultiLabelBinarizer(classes=all_labels_sorted)
    mlb.fit(df["true_set"])
    label_to_class = {lbl: lbl.split(".")[0] for lbl in all_labels_sorted}

    results: dict = {}
    for col in prediction_cols:
        y_true = mlb.transform(df["true_set"])
        y_pred = mlb.transform(df[col + "_pred"])
        support = y_true.sum(axis=0)

        mcc_list = [matthews_corrcoef(y_true[:, i], y_pred[:, i]) for i in range(len(all_labels_sorted))]
        prec_list = precision_score(y_true, y_pred, average=None, zero_division=0)
        rec_list = recall_score(y_true, y_pred, average=None, zero_division=0)

        overall_mcc = float(np.average(mcc_list, weights=support))
        overall_prec = float(np.average(prec_list, weights=support))
        overall_rec = float(np.average(rec_list, weights=support))
        coverage = float((df[col + "_pred"].apply(len) > 0).mean())

        class_mcc: dict[str, float] = {}
        for cls in sorted(set(label_to_class.values())):
            idx = [i for i, lbl in enumerate(all_labels_sorted) if label_to_class[lbl] == cls]
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


def find_gemma_tsv(gemma_split_dir: str) -> str:
    """Return the first TSV matching the Gemma4EC predictions pattern."""
    import glob
    pattern = f"{gemma_split_dir}/results/predictions_pred.EC_NUMBER.*.tsv"
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f"No Gemma4EC TSV found matching: {pattern}")
    return matches[0]


def main() -> None:
    args = parse_args()

    splits_cfg = {
        "Scaffold": {
            "test_tsv": f"{args.rhea_data_dir}/Scaffold/test.tsv",
            "becpred":  f"{args.becpred_dir}/Scaffold_becpred.csv",
            "theia":    f"{args.theia_dir}/scaffold_ec123_test_predictions.csv",
            "gemma_dir": f"{args.gemma_dir}/Rhea_Scaffold",
        },
        "Stratified": {
            "test_tsv": f"{args.rhea_data_dir}/Stratified/test.tsv",
            "becpred":  f"{args.becpred_dir}/Stratified_becpred.csv",
            "theia":    f"{args.theia_dir}/stratified_ec123_test_predictions.csv",
            "gemma_dir": f"{args.gemma_dir}/Rhea_Stratified",
        },
        "Time": {
            "test_tsv": f"{args.rhea_data_dir}/Time/test.tsv",
            "becpred":  f"{args.becpred_dir}/Time_becpred.csv",
            "theia":    f"{args.theia_dir}/time_ec123_test_predictions.csv",
            "gemma_dir": f"{args.gemma_dir}/Rhea_Time",
        },
    }

    all_records: list[dict] = []

    for split, paths in splits_cfg.items():
        print(f"\n=== {split} ===")

        test = pd.read_csv(paths["test_tsv"], sep="\t")
        test.columns = [c.strip() for c in test.columns]
        test = test.rename(columns={
            "REACTION_ID": "reaction_id",
            "REACTION_SMILES": "smiles",
            "EC_NUMBER": "ec_number",
        })
        test["true_set"] = test["ec_number"].apply(parse_ecs)
        test = test[test["true_set"].apply(len) > 0].copy()

        bec = pd.read_csv(paths["becpred"])
        bec.columns = [c.strip() for c in bec.columns]
        bec = bec.rename(columns={"Reaction ID": "reaction_id", "Prediction": "BEC-Pred"})
        bec["reaction_id"] = bec["reaction_id"].astype(str)
        test["reaction_id"] = test["reaction_id"].astype(str)
        df = test.merge(bec[["reaction_id", "BEC-Pred"]], on="reaction_id", how="left")

        theia = pd.read_csv(paths["theia"])
        theia = theia[["rxn_smiles", "top1_ec"]].rename(
            columns={"rxn_smiles": "smiles", "top1_ec": "Theia"}
        )
        df = df.merge(theia, on="smiles", how="left")

        gemma_tsv = find_gemma_tsv(paths["gemma_dir"])
        gemma = pd.read_csv(gemma_tsv, sep="\t")
        gemma = gemma[["smiles", "pred_fragment"]].rename(columns={"pred_fragment": "Gemma4EC"})
        df = df.merge(gemma, on="smiles", how="left")

        for col in ["BEC-Pred", "Theia", "Gemma4EC"]:
            df[col] = df[col].fillna("")
            df[col + "_pred"] = df[col].apply(parse_ecs)

        n = len(df)
        for col in ["BEC-Pred", "Theia", "Gemma4EC"]:
            covered = (df[col + "_pred"].apply(len) > 0).sum()
            print(f"  {col}: {covered}/{n} reactions with predictions")

        results = evaluate_methods(df, ["BEC-Pred", "Theia", "Gemma4EC"])

        for method, res in results.items():
            print(
                f"  {method}: MCC={res['overall_mcc']:.4f}, "
                f"Prec={res['overall_precision']:.4f}, "
                f"Rec={res['overall_recall']:.4f}, "
                f"Cov={res['coverage']:.4f}"
            )
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
    out.to_csv(args.output, index=False)
    print(f"\nSaved to {args.output}")

    print("\n=== OVERALL MCC SUMMARY ===")
    summary = out.drop_duplicates(subset=["split", "method"])[
        ["split", "method", "overall_mcc", "coverage", "total_support"]
    ]
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
