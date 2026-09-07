"""
Stratified random 90/10 train/test split by EC subsubclass.

Each EC subsubclass is represented in both splits at the same 90/10
ratio, ensuring balanced class coverage.

Splitting happens at the REACTION_ID level, not the row level: a reaction
with multiple EC rows (broad-specificity enzymes, one row per EC) is kept
entirely in train or entirely in test, using its first EC row as the
stratification label. Splitting per-row instead would let the same
reaction's SMILES leak across train and test under different EC labels.
"""

import argparse
import os

import pandas as pd
from sklearn.model_selection import train_test_split


def get_ec_subsubclass(ec_number: str) -> str | None:
    try:
        parts = str(ec_number).split('.')
        if len(parts) >= 3:
            return '.'.join(parts[:3])
        return str(ec_number)
    except Exception:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stratified random 90/10 train/test split of Rhea reaction data."
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to input TSV with columns REACTION_ID, REACTION_SMILES, EC_NUMBER."
    )
    parser.add_argument(
        "--output_dir", default=".",
        help="Directory to write train.tsv and test.tsv (default: current directory)."
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed (default: 42)."
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading data...")
    df = pd.read_csv(args.input, sep='\t')
    print(f"Total reactions: {len(df)}")

    print("\nExtracting EC subsubclasses...")
    df['ec_subsubclass'] = df['EC_NUMBER'].apply(get_ec_subsubclass)
    df_clean = df.dropna(subset=['ec_subsubclass'])
    print(f"Reactions with valid EC: {len(df_clean)}")

    n_multi_ec_ids = (df_clean.groupby('REACTION_ID')['EC_NUMBER'].nunique() > 1).sum()
    if n_multi_ec_ids:
        print(f"{n_multi_ec_ids} REACTION_ID(s) carry more than one EC row -- "
              f"splitting by REACTION_ID group so these stay entirely in one split.")

    # One representative row per REACTION_ID drives the group-level stratified
    # split; all EC rows for that REACTION_ID then follow the same assignment.
    groups = df_clean.drop_duplicates(subset=['REACTION_ID'], keep='first')

    print("\nPerforming stratified random split (90% train, 10% test)...")
    print("-" * 80)

    class_counts = groups['ec_subsubclass'].value_counts()
    singleton_classes = class_counts[class_counts < 2].index
    groups_singleton = groups[groups['ec_subsubclass'].isin(singleton_classes)]
    groups_splittable = groups[~groups['ec_subsubclass'].isin(singleton_classes)]
    if len(groups_singleton) > 0:
        print(
            f"Note: {len(singleton_classes)} EC subsubclass(es) have only 1 reaction "
            f"({len(groups_singleton)} reactions total) — cannot be split, routed entirely to train."
        )

    groups_train, groups_test = train_test_split(
        groups_splittable,
        test_size=0.1,
        stratify=groups_splittable['ec_subsubclass'],
        random_state=args.seed,
    )
    train_ids = set(groups_train['REACTION_ID']) | set(groups_singleton['REACTION_ID'])
    test_ids = set(groups_test['REACTION_ID'])

    train_df = df_clean[df_clean['REACTION_ID'].isin(train_ids)]
    test_df = df_clean[df_clean['REACTION_ID'].isin(test_ids)]

    for ec_class, group in df_clean.groupby('ec_subsubclass'):
        n_total = len(group)
        n_train = len(train_df[train_df['ec_subsubclass'] == ec_class])
        n_test = len(test_df[test_df['ec_subsubclass'] == ec_class])
        print(
            f"EC {ec_class}: total {n_total:5d} → "
            f"train: {n_train:5d} ({n_train/n_total*100:.1f}%), "
            f"test: {n_test:5d} ({n_test/n_total*100:.1f}%)"
        )

    print("\n" + "=" * 80)
    print(f"TRAIN: {len(train_df):6d} reactions ({len(train_df)/len(df_clean)*100:.1f}%)")
    print(f"TEST:  {len(test_df):6d} reactions ({len(test_df)/len(df_clean)*100:.1f}%)")
    print("=" * 80)

    cols = ['REACTION_ID', 'REACTION_SMILES', 'EC_NUMBER']
    train_path = os.path.join(args.output_dir, 'train.tsv')
    test_path = os.path.join(args.output_dir, 'test.tsv')
    train_df[cols].to_csv(train_path, sep='\t', index=False)
    test_df[cols].to_csv(test_path, sep='\t', index=False)

    print(f"\nFiles saved: {train_path}, {test_path}")


if __name__ == "__main__":
    main()
