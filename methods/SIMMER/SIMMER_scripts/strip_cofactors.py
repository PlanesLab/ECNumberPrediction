"""
Strip cofactor components from a reaction SMILES column and write a
no-cofactor variant of the input file.

Reactions where either side becomes empty after stripping (i.e. the whole
side was made of cofactors, such as a bare metal-ion redox reaction) are
dropped and logged separately, since there is nothing left to fingerprint.

Usage:
    python strip_cofactors.py \
        --input data/Splits-Rhea/Stratified/train.tsv \
        --output data/Splits-Rhea/Stratified/train_nocofactor.tsv \
        --smiles-col REACTION_SMILES \
        --sep "\t"
"""

import argparse
import pandas as pd

from cofactors import strip_cofactor_components


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--smiles-col", default="REACTION_SMILES")
    parser.add_argument("--sep", default="\t")
    parser.add_argument("--dropped-log", default=None,
                         help="Optional path for a CSV log of removed cofactor components")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep=args.sep)

    kept_rows = []
    dropped_rows = []
    log_rows = []

    for _, row in df.iterrows():
        left, right = row[args.smiles_col].split(">>")
        left_kept, left_dropped = strip_cofactor_components(left.strip())
        right_kept, right_dropped = strip_cofactor_components(right.strip())

        for comp in left_dropped:
            log_rows.append({"reaction": row.get(df.columns[0]), "side": "left", "cofactor_smiles": comp})
        for comp in right_dropped:
            log_rows.append({"reaction": row.get(df.columns[0]), "side": "right", "cofactor_smiles": comp})

        if not left_kept or not right_kept:
            dropped_rows.append(row)
            continue

        new_row = row.copy()
        new_row[args.smiles_col] = f"{left_kept}>>{right_kept}"
        kept_rows.append(new_row)

    kept_df = pd.DataFrame(kept_rows)
    kept_df.to_csv(args.output, sep=args.sep, index=False)

    dropped_log_path = args.dropped_log or (args.output + ".dropped_reactions.csv")
    pd.DataFrame(dropped_rows).to_csv(dropped_log_path, index=False)

    cofactor_log_path = args.output + ".cofactor_components.csv"
    pd.DataFrame(log_rows).to_csv(cofactor_log_path, index=False)

    print(f"Input reactions:            {len(df)}")
    print(f"Kept (no-cofactor) reactions: {len(kept_df)} -> {args.output}")
    print(f"Dropped (both/one side empty after stripping): {len(dropped_rows)} -> {dropped_log_path}")
    print(f"Cofactor components removed log: {len(log_rows)} -> {cofactor_log_path}")


if __name__ == "__main__":
    main()
