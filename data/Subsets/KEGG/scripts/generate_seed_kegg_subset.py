"""
Generates a seeded, EC-stratified resample of the KEGG-1.8K Case 1 test set.

Pool = data/Splits-DBs/KEGG/{train,test}.tsv combined (7782 reactions). Draws
a stratified subset the same size as the original test set (1866 reactions)
for a given --seed, using stratified_split.py's EC-subsubclass pattern.

Writes two row-order-aligned outputs so any Case 1 run script can be
repointed at a seed by swapping these two paths:
  <output_dir>/kegg_reactions_current_test.csv   (drop-in replacement for the test CSV)
  <output_dir>/8:2KEGGTest_canonicalized.txt      (one canonical SMILES per line, same order)
"""

import argparse
import os
import sys

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from sklearn.model_selection import train_test_split

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'Preprocessing'))
from canonicalize_rxn_SMILES import canonicalize_reaction_smiles  # noqa: E402

POOL_TRAIN = "data/Splits-DBs/KEGG/train.tsv"
POOL_TEST = "data/Splits-DBs/KEGG/test.tsv"
ORIGINAL_TEST_SIZE = 1866


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seeded stratified resample of the KEGG-1.8K test set.")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--pool_train", default=POOL_TRAIN)
    parser.add_argument("--pool_test", default=POOL_TEST)
    parser.add_argument("--sample_size", type=int, default=ORIGINAL_TEST_SIZE)
    parser.add_argument("--output_dir", required=True)
    return parser.parse_args()


def get_ec_subsubclass(ec_number: str):
    parts = str(ec_number).split('.')
    return '.'.join(parts[:3]) if len(parts) >= 3 else str(ec_number)


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    pool = pd.concat([
        pd.read_csv(args.pool_train, sep='\t', dtype=str),
        pd.read_csv(args.pool_test, sep='\t', dtype=str),
    ], ignore_index=True)
    pool = pool.rename(columns={'ec': 'EC Number'})
    print(f"Pool size: {len(pool)}")

    pool['ec_subsubclass'] = pool['EC Number'].apply(get_ec_subsubclass)
    class_counts = pool['ec_subsubclass'].value_counts()
    singleton_classes = class_counts[class_counts < 2].index
    df_singleton = pool[pool['ec_subsubclass'].isin(singleton_classes)]
    df_splittable = pool[~pool['ec_subsubclass'].isin(singleton_classes)]

    frac = args.sample_size / len(pool)
    _, sample_df = train_test_split(
        df_splittable,
        test_size=frac,
        stratify=df_splittable['ec_subsubclass'],
        random_state=args.seed,
    )
    # Singleton EC classes can't be stratified; include proportionally at random.
    if len(df_singleton) > 0:
        n_singleton_sample = max(0, round(len(df_singleton) * frac))
        sample_df = pd.concat([
            sample_df,
            df_singleton.sample(n=min(n_singleton_sample, len(df_singleton)), random_state=args.seed),
        ])

    sample_df = sample_df.sample(frac=1, random_state=args.seed).reset_index(drop=True)
    sample_df = sample_df.drop(columns=['ec_subsubclass'])
    print(f"Seed {args.seed}: sampled {len(sample_df)} reactions")

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
