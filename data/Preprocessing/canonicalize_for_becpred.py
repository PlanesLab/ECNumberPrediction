"""
Canonicalizes a Rhea seed split's reaction SMILES for BEC-Pred.

Falls back to the original SMILES on RDKit failure instead of dropping the row.
BEC-Pred's queries.txt is positional, so a dropped row desyncs predictions from
ground truth.

Writes train.tsv/test.tsv (reaction_smiles + rxn columns canonicalized) and
queries.txt, regenerated from the canonicalized test.tsv.
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from canonicalize_rxn_SMILES import canonicalize_reaction_smiles


def canonicalize_with_fallback(smiles: str) -> str:
    canonical = canonicalize_reaction_smiles(smiles)
    return canonical if canonical is not None else smiles


def main() -> None:
    parser = argparse.ArgumentParser(description="Canonicalize a prepared Rhea split's SMILES for BEC-Pred.")
    parser.add_argument("--train", required=True, help="prepared/train.tsv")
    parser.add_argument("--test", required=True, help="prepared/test.tsv")
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    n_total = 0
    n_fallback = 0

    def process(path: str, out_name: str) -> pd.DataFrame:
        nonlocal n_total, n_fallback
        df = pd.read_csv(path, sep="\t", dtype=str)
        n_total += len(df)
        canonical = df["reaction_smiles"].apply(canonicalize_with_fallback)
        n_fallback += df["reaction_smiles"].apply(lambda s: canonicalize_reaction_smiles(s) is None).sum()
        df["reaction_smiles"] = canonical
        if "rxn" in df.columns:
            df["rxn"] = canonical
        out_path = os.path.join(args.output_dir, out_name)
        df.to_csv(out_path, sep="\t", index=False)
        return df

    process(args.train, "train.tsv")
    test_df = process(args.test, "test.tsv")

    queries_path = os.path.join(args.output_dir, "queries.txt")
    with open(queries_path, "w") as f:
        for s in test_df["reaction_smiles"]:
            f.write(str(s) + "\n")

    print(f"Canonicalized {n_total} rows ({n_fallback} could not be canonicalized, kept original SMILES)")
    print(f"Wrote '{args.output_dir}/train.tsv', 'test.tsv', 'queries.txt'")


if __name__ == "__main__":
    main()
