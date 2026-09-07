"""
Generate a seeded, EC-stratified 90/10 train/test split of the MetaNetX
Case 2 reaction pool.

data/Splits-DBs/MetaNetX/{train,test}.tsv were pure carry-overs from before
this repo existed (no generating script was ever committed) at a fixed
90/10 ratio (34047/3784 -- confirmed by direct measurement, despite the
README's "80/20" claim). This script combines them back into one pool and
draws a fresh stratified 90/10 split per --seed, reusing the same
stratify-by-EC-subsubclass pattern as data/Splits-Rhea/scripts/stratified_split.py
and data/Subsets/KEGG/scripts/generate_seed_kegg_subset.py.

NOTE: the original train.tsv and test.tsv disagree on the SMILES column
name ('rxn' vs 'reaction_smiles') -- a real inconsistency in the existing
data. This script standardizes on 'reaction_smiles' for both outputs (the
name run_SIMMER.sh's Case 2 block already expects for both its train.tsv
and test.tsv reads), so no downstream script needs further column-name
patching for that specific issue.

Writes, per seed:
  <output_dir>/train.tsv    (90%, columns reaction_id, reaction_smiles, substrates_products, ec)
  <output_dir>/test.tsv     (10%, same columns)
  <output_dir>/queries.txt  (test.tsv's reaction_smiles, one per line, same row order as test.tsv --
                              required by Theia/BEC-Pred/CLAIRE's positional query scripts)
"""

import argparse
import os

import pandas as pd
from sklearn.model_selection import train_test_split

POOL_TRAIN = "data/Splits-DBs/MetaNetX/train.tsv"
POOL_TEST = "data/Splits-DBs/MetaNetX/test.tsv"
TEST_FRACTION = 0.10  # matches the actual measured ratio of the existing baseline split


def get_ec_subsubclass(ec_number: str):
    parts = str(ec_number).split('.')
    return '.'.join(parts[:3]) if len(parts) >= 3 else str(ec_number)


def load_pool(pool_train: str, pool_test: str) -> pd.DataFrame:
    train = pd.read_csv(pool_train, sep='\t', dtype=str).rename(columns={'rxn': 'reaction_smiles'})
    test = pd.read_csv(pool_test, sep='\t', dtype=str).rename(columns={'rxn': 'reaction_smiles'})
    pool = pd.concat([train, test], ignore_index=True)
    return pool[['reaction_id', 'reaction_smiles', 'substrates_products', 'ec']]


def main() -> None:
    parser = argparse.ArgumentParser(description="Seeded 90/10 stratified split of the MetaNetX pool.")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--pool_train", default=POOL_TRAIN)
    parser.add_argument("--pool_test", default=POOL_TEST)
    parser.add_argument("--test_fraction", type=float, default=TEST_FRACTION)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    pool = load_pool(args.pool_train, args.pool_test)
    print(f"Pool size: {len(pool)}")

    pool['ec_subsubclass'] = pool['ec'].apply(get_ec_subsubclass)
    class_counts = pool['ec_subsubclass'].value_counts()
    singleton_classes = class_counts[class_counts < 2].index
    df_singleton = pool[pool['ec_subsubclass'].isin(singleton_classes)]
    df_splittable = pool[~pool['ec_subsubclass'].isin(singleton_classes)]
    if len(df_singleton) > 0:
        print(f"{len(singleton_classes)} singleton EC subsubclass(es) ({len(df_singleton)} rows) routed to train")

    train_df, test_df = train_test_split(
        df_splittable,
        test_size=args.test_fraction,
        stratify=df_splittable['ec_subsubclass'],
        random_state=args.seed,
    )
    train_df = pd.concat([train_df, df_singleton]).drop(columns=['ec_subsubclass'])
    test_df = test_df.drop(columns=['ec_subsubclass']).sample(frac=1, random_state=args.seed).reset_index(drop=True)

    print(f"Seed {args.seed}: train={len(train_df)} ({len(train_df)/len(pool)*100:.1f}%), "
          f"test={len(test_df)} ({len(test_df)/len(pool)*100:.1f}%)")

    train_path = os.path.join(args.output_dir, "train.tsv")
    test_path = os.path.join(args.output_dir, "test.tsv")
    queries_path = os.path.join(args.output_dir, "queries.txt")

    # Some consumers (e.g. theia's encode_split_data.py) still expect the
    # pre-reorg 'rxn' column name instead of 'reaction_smiles' -- write both
    # rather than chasing every consumer's own assumption.
    train_df.insert(2, 'rxn', train_df['reaction_smiles'])
    test_df.insert(2, 'rxn', test_df['reaction_smiles'])

    train_df.to_csv(train_path, sep='\t', index=False)
    test_df.to_csv(test_path, sep='\t', index=False)
    with open(queries_path, "w") as f:
        for s in test_df['reaction_smiles']:
            f.write(str(s) + "\n")

    print(f"Wrote '{train_path}', '{test_path}', '{queries_path}' (queries.txt row-order aligned with test.tsv)")


if __name__ == "__main__":
    main()
