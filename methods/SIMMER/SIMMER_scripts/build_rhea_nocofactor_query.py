"""
Build a SIMMER2 query CSV (5 columns: reaction, left_comp, right_comp,
left_smiles, right_smiles) from a no-cofactor Rhea test split, plus a
ground-truth CSV mapping each query id to its true EC number for later
evaluation.

Query ids are "{REACTION_ID}_{EC_NUMBER}_Q" -- guaranteed disjoint from the
"{REACTION_ID}_{EC_NUMBER}" ids used in the training DB (see
build_rhea_nocofactor_db.py), since Rhea REACTION_IDs are not unique and can
otherwise collide between the train and test splits.
"""

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-query", required=True)
    parser.add_argument("--output-truth", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    df["EC_NUMBER"] = df["EC_NUMBER"].astype(str).str.strip()
    df = df[df["EC_NUMBER"].apply(lambda x: x and x.lower() != "nan")].reset_index(drop=True)

    # Hyphens (not underscores) keep the id intact through ec_predictions.py's
    # `filename.split('_')[0]` reaction-name extraction downstream, and the
    # "-Q" suffix guarantees no collision with build_rhea_nocofactor_db.py ids.
    query_id = df["REACTION_ID"].astype(str) + "-" + df["EC_NUMBER"] + "-Q"
    assert query_id.is_unique, "query ids are not unique"

    left, right = zip(*df["REACTION_SMILES"].str.split(">>"))
    out = pd.DataFrame({
        "reaction": query_id,
        "left_comp": [s.strip() for s in left],
        "right_comp": [s.strip() for s in right],
        "left_smiles": [s.strip() for s in left],
        "right_smiles": [s.strip() for s in right],
    })
    out.to_csv(args.output_query, index=False)

    truth = pd.DataFrame({
        "reaction": query_id,
        "REACTION_ID": df["REACTION_ID"],
        "true_EC": df["EC_NUMBER"],
    })
    truth.to_csv(args.output_truth, index=False)

    print(f"Query reactions: {len(out)} -> {args.output_query}")
    print(f"Ground truth:    {len(truth)} -> {args.output_truth}")


if __name__ == "__main__":
    main()
