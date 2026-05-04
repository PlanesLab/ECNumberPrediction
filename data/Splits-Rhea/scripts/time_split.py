"""
Time-based train/test split using two Rhea database snapshots.

Reactions present in the older snapshot become the training set.
Reactions added in the newer snapshot become the test set. This
simulates a real-world temporal evaluation where models are trained
on historical data and tested on future reactions.
"""

import argparse
import os

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Time-based train/test split: reactions from an older DB snapshot "
            "are used for training; new reactions form the test set."
        )
    )
    parser.add_argument(
        "--old_db", required=True,
        help="Path to the older database TSV (e.g., Rhea 2018 snapshot)."
    )
    parser.add_argument(
        "--new_db", required=True,
        help="Path to the current (full) database TSV."
    )
    parser.add_argument(
        "--output_dir", default=".",
        help="Directory to write train.tsv and test.tsv (default: current directory)."
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading old database...")
    df_old = pd.read_csv(args.old_db, sep='\t')
    print(f"Old DB reactions: {len(df_old)}")

    print("\nLoading new database...")
    df_new = pd.read_csv(args.new_db, sep='\t')
    print(f"New DB reactions: {len(df_new)}")

    old_ids = set(df_old['REACTION_ID'])
    print(f"\nUnique REACTION_IDs in old DB: {len(old_ids)}")

    train_df = df_new[df_new['REACTION_ID'].isin(old_ids)]
    test_df = df_new[~df_new['REACTION_ID'].isin(old_ids)]

    print("\n" + "=" * 80)
    print("TEMPORAL SPLIT RESULTS:")
    print("-" * 80)
    print(f"TRAIN (old reactions): {len(train_df):6d} ({len(train_df)/len(df_new)*100:.1f}%)")
    print(f"TEST  (new reactions): {len(test_df):6d} ({len(test_df)/len(df_new)*100:.1f}%)")
    print("=" * 80)

    train_path = os.path.join(args.output_dir, 'train.tsv')
    test_path = os.path.join(args.output_dir, 'test.tsv')
    train_df.to_csv(train_path, sep='\t', index=False)
    test_df.to_csv(test_path, sep='\t', index=False)

    missed = len(old_ids - set(train_df['REACTION_ID']))
    print(f"\nFiles saved: {train_path}, {test_path}")
    print(f"Reactions in old DB not found in new DB: {missed}")


if __name__ == "__main__":
    main()
