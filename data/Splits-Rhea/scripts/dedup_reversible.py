"""
Collapse Rhea's reversible reaction pairs (forward/reverse SMILES of the
same underlying reaction, e.g. RHEA LR/RL entries) into a single row each,
so no reaction is represented twice under the same EC number.

For every reaction the reverse SMILES (products>>reactants) is computed and
paired against the rest of the dataset. If the exact reverse (same EC) is
also present, only one of the two directions is kept (the lower
REACTION_ID, for reproducibility). If no reverse partner exists, the row is
kept as-is. Rows sharing a REACTION_ID with a *different* EC (multi-label
reactions) are left untouched — only forward/reverse duplication is
collapsed.
"""

import argparse
from typing import Optional

import pandas as pd


def reverse_smiles(reaction_smiles: str) -> Optional[str]:
    parts = reaction_smiles.split(">>")
    if len(parts) != 2:
        return None
    lhs, rhs = parts
    return f"{rhs}>>{lhs}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Deduplicate reversible Rhea reaction pairs (keep one direction per pair)."
    )
    parser.add_argument("--input", required=True, help="Path to input TSV with REACTION_ID, REACTION_SMILES, EC column.")
    parser.add_argument("--output", required=True, help="Path to write the deduplicated TSV.")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    smiles_col, ec_col = df.columns[1], df.columns[2]
    print(f"Input rows: {len(df)}")

    df = df.drop_duplicates()
    print(f"After exact-row dedup: {len(df)}")

    df["_rev_smiles"] = df[smiles_col].apply(reverse_smiles)
    df["_pair_key"] = df.apply(
        lambda r: tuple(sorted([r[smiles_col], r["_rev_smiles"]])), axis=1
    )

    before = len(df)
    df = (
        df.sort_values("REACTION_ID")
        .drop_duplicates(subset=["_pair_key", ec_col], keep="first")
        .drop(columns=["_rev_smiles", "_pair_key"])
        .sort_values("REACTION_ID")
        .reset_index(drop=True)
    )
    print(f"After reversible-pair dedup: {len(df)} (removed {before - len(df)})")

    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
