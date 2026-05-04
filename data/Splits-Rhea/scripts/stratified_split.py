"""
Stratified random 90/10 train/test split by EC subsubclass.

Each EC subsubclass is represented in both splits at the same 90/10
ratio, ensuring balanced class coverage.
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

    print("\nPerforming stratified random split (90% train, 10% test)...")
    print("-" * 80)

    train_df, test_df = train_test_split(
        df_clean,
        test_size=0.1,
        stratify=df_clean['ec_subsubclass'],
        random_state=args.seed,
    )

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
