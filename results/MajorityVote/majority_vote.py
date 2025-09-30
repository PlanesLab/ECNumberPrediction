import argparse
import pandas as pd
from collections import Counter

# === Functions ===
def collapse_to_third_level(predictions):
    if pd.isna(predictions) or predictions == "" or str(predictions).lower() == "nan":
        return None
    preds = str(predictions).split(";")
    collapsed = []
    for p in preds:
        sub_preds = [x.strip() for x in p.split("|") if x.strip() != ""]
        sub_preds = list({ ".".join(x.split(".")[:3]) for x in sub_preds })
        collapsed.append("|".join(sub_preds))
    return ";".join(collapsed)

def extract_top1(prediction):
    if not prediction:
        return []
    preds = prediction.split(";")
    first = preds[0]
    return [x.strip() for x in first.split("|") if x.strip() != ""]

def extract_top5_rstyle(prediction):
    """Extract top5 like R: pipe = same level, semicolon = next level, weight per group"""
    if not prediction:
        return []
    groups = prediction.split(";")
    top5 = []
    for grp in groups[:5]:  # limit to first 5 semicolon groups like R
        alternatives = [x.strip() for x in grp.split("|") if x.strip() != ""]
        top5.extend(alternatives)
    return top5[:5]  # ensure max 5 predictions

def majority_vote_rstyle(preds_dict, mode="top1"):
    weighted_votes = Counter()
    for method, preds in preds_dict.items():
        if not preds:
            continue
        if mode == "top1":
            for p in preds:
                weighted_votes[p] += 5
        elif mode == "top5":
            for i, p in enumerate(preds):
                # weight = max(6 - group_index, 1) → emulate R
                weighted_votes[p] += max(6 - i, 1)
    if not weighted_votes:
        return None
    return weighted_votes.most_common(1)[0][0]

# === Main ===
parser = argparse.ArgumentParser(description="Compute Top1/Top5 majority EC predictions")
parser.add_argument("--input_csv", type=str, required=True, help="Path to merged_output CSV")
parser.add_argument("--drug", type=str, required=True, help="Drug name or 'all'")
parser.add_argument("--output_csv", type=str, default=None, help="Optional output CSV for all drugs")
args = parser.parse_args()

# Load CSV
df = pd.read_csv(args.input_csv)
df = df.fillna("nan")

# Columns containing EC predictions
ec_methods = ["E-zyme1","E-zyme2", "BridgIT", "SelenzymeRF","SIMMER", "Theia_ECREACT" ,"BEC-Pred","CLAIRE"]

# Collapse predictions to third level
for m in ec_methods:
    df[m] = df[m].apply(collapse_to_third_level)

# Methods used for majority vote
majority_methods = ["Theia_ECREACT","BEC-Pred","SIMMER","SelenzymeRF"]

# Function to compute majority votes for a row
def compute_majority_rstyle(row):
    top1_preds_dict = {m: extract_top1(row[m]) for m in majority_methods}
    top5_preds_dict = {m: extract_top5_rstyle(row[m]) for m in majority_methods}
    return pd.Series({
        "majority_top1": majority_vote_rstyle(top1_preds_dict, mode="top1"),
        "majority_top5": majority_vote_rstyle(top5_preds_dict, mode="top5")
    })

# === Filter and compute ===
if args.drug.lower() != "all":
    drug_row = df[df['drug'] == args.drug]
    if drug_row.empty:
        print(f"No data found for drug: {args.drug}")
        exit()
    majority_results = compute_majority_rstyle(drug_row.iloc[0])
    print(f"Drug: {args.drug}")
    print(f"Top 1 majority vote EC: {majority_results['majority_top1']}")
    print(f"Top 5 majority vote EC: {majority_results['majority_top5']}")
else:
    # All drugs
    majority_df = df.copy()
    majority_df[["majority_top1","majority_top5"]] = majority_df.apply(compute_majority_rstyle, axis=1)
    print(majority_df[["drug","majority_top1","majority_top5"]])
    if args.output_csv:
        majority_df[["drug","majority_top1","majority_top5"]].to_csv(args.output_csv, index=False)
        print(f"Saved majority votes for all drugs to {args.output_csv}")