"""
Generates a seeded bootstrap sample (default n=100, with replacement) of the
full KEGG reaction pool, for Case 1's bootstrap re-run.

Pool = data/Splits-DBs/KEGG/{train,test}.tsv combined (7782 reactions), same
pool generate_seed_kegg_subset.py draws from -- but sampled uniformly with
replacement rather than stratified without replacement.

Writes the same two row-order-aligned outputs as generate_seed_kegg_subset.py:
  <output_dir>/kegg_reactions_current_test.csv   (drop-in replacement for the test CSV)
  <output_dir>/8:2KEGGTest_canonicalized.txt      (one canonical SMILES per line, same order)
"""

import argparse
import os
import sys

import pandas as pd
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'Preprocessing'))
from canonicalize_rxn_SMILES import canonicalize_reaction_smiles  # noqa: E402

POOL_TRAIN = "data/Splits-DBs/KEGG/train.tsv"
POOL_TEST = "data/Splits-DBs/KEGG/test.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seeded bootstrap (with-replacement) sample of the full KEGG pool.")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--pool_train", default=POOL_TRAIN)
    parser.add_argument("--pool_test", default=POOL_TEST)
    parser.add_argument("--sample_size", type=int, default=100)
    parser.add_argument("--output_dir", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    pool = pd.concat([
        pd.read_csv(args.pool_train, sep='\t', dtype=str),
        pd.read_csv(args.pool_test, sep='\t', dtype=str),
    ], ignore_index=True)
    pool = pool.rename(columns={'ec': 'EC Number'})
    print(f"Pool size: {len(pool)}")

    sample_df = pool.sample(n=args.sample_size, replace=True, random_state=args.seed).reset_index(drop=True)
    n_dup = len(sample_df) - sample_df['Reaction ID'].nunique()
    print(f"Seed {args.seed}: bootstrap-sampled {len(sample_df)} reactions with replacement "
          f"({n_dup} duplicate draws)")

    canon_smiles = sample_df['Isomeric_SMILES'].astype(str).apply(canonicalize_reaction_smiles)
    n_failed = canon_smiles.isna().sum()
    if n_failed > 0:
        print(f"WARNING: {n_failed} reactions failed canonicalization, dropping to keep row alignment")
        keep_mask = canon_smiles.notna()
        sample_df = sample_df[keep_mask].reset_index(drop=True)
        canon_smiles = canon_smiles[keep_mask].reset_index(drop=True)

    csv_path = os.path.join(args.output_dir, "kegg_reactions_current_test.csv")
    txt_path = os.path.join(args.output_dir, "8:2KEGGTest_canonicalized.txt")
    sample_df.to_csv(csv_path, index=False)
    with open(txt_path, "w") as f:
        for s in canon_smiles:
            f.write(s + "\n")

    print(f"Wrote {len(sample_df)} rows to '{csv_path}' and '{txt_path}' (row-order aligned)")


if __name__ == "__main__":
    main()
