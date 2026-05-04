"""
Scaffold-aware 90/10 train/test split using Bemis-Murcko scaffolds.

Splits are stratified per EC subsubclass so that scaffolds in train and
test do not overlap, reducing data leakage from structurally similar
molecules.
"""

import argparse
import os

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def get_bemis_murcko_scaffold(smiles: str) -> str | None:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except Exception:
        return None


def extract_product_smiles(reaction_smiles: str) -> str | None:
    try:
        parts = reaction_smiles.split('>>')
        if len(parts) == 2:
            return parts[1].split('.')[0].strip()
        return None
    except Exception:
        return None


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
        description="Scaffold-aware 90/10 train/test split of Rhea reaction data."
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

    print("\nExtracting Bemis-Murcko scaffolds...")
    df['product_smiles'] = df['REACTION_SMILES'].apply(extract_product_smiles)
    df['scaffold'] = df['product_smiles'].apply(get_bemis_murcko_scaffold)
    df['ec_subsubclass'] = df['EC_NUMBER'].apply(get_ec_subsubclass)

    df_clean = df.dropna(subset=['scaffold', 'ec_subsubclass'])
    print(f"Reactions with valid scaffolds: {len(df_clean)}")

    print("\nPerforming scaffold-aware split (90% train, 10% test)...")
    print("-" * 80)

    rng = np.random.default_rng(args.seed)
    train_indices: list[int] = []
    test_indices: list[int] = []

    for ec_class, group in df_clean.groupby('ec_subsubclass'):
        scaffolds = group['scaffold'].unique()
        n_train = max(1, int(0.9 * len(scaffolds)))

        shuffled = rng.permutation(scaffolds)
        train_scaffolds = set(shuffled[:n_train])
        test_scaffolds = set(shuffled[n_train:])

        group_train = group[group['scaffold'].isin(train_scaffolds)]
        group_test = group[group['scaffold'].isin(test_scaffolds)]

        train_indices.extend(group_train.index.tolist())
        test_indices.extend(group_test.index.tolist())

        print(
            f"EC {ec_class}: {len(scaffolds):4d} scaffolds → "
            f"train: {len(train_scaffolds):4d} ({len(group_train):5d} rxns), "
            f"test: {len(test_scaffolds):4d} ({len(group_test):5d} rxns)"
        )

    train_df = df_clean.loc[train_indices]
    test_df = df_clean.loc[test_indices]

    print("\n" + "=" * 80)
    print(f"TRAIN: {len(train_df):6d} reactions ({len(train_df)/len(df_clean)*100:.1f}%)")
    print(f"TEST:  {len(test_df):6d} reactions ({len(test_df)/len(df_clean)*100:.1f}%)")
    print("=" * 80)

    cols = ['REACTION_ID', 'REACTION_SMILES', 'EC_NUMBER']
    train_path = os.path.join(args.output_dir, 'train.tsv')
    test_path = os.path.join(args.output_dir, 'test.tsv')
    train_df[cols].to_csv(train_path, sep='\t', index=False)
    test_df[cols].to_csv(test_path, sep='\t', index=False)

    train_scf = set(train_df['scaffold'])
    test_scf = set(test_df['scaffold'])
    print(f"\nUnique scaffolds in train: {len(train_scf)}")
    print(f"Unique scaffolds in test:  {len(test_scf)}")
    print(f"Scaffold overlap:          {len(train_scf & test_scf)}")
    print(f"\nFiles saved: {train_path}, {test_path}")


if __name__ == "__main__":
    main()
