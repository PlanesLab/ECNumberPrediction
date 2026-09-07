"""
Renames a Rhea seed split file (uppercase REACTION_ID/REACTION_SMILES/EC_NUMBER)
to the lowercase MetaNetX-style schema (reaction_id, reaction_smiles, rxn, ec)
that BEC-Pred/Theia/CLAIRE hardcode column names for. Optionally emits a
queries.txt alongside it.

SIMMER takes column names as CLI flags and reads the raw uppercase files
directly, so it doesn't need this script.

Writes:
  <output>            train.tsv or test.tsv, columns reaction_id, reaction_smiles, rxn, ec
  <queries_output>    (optional) reaction_smiles, one per line, same row order. Only
                       meaningful when --input is a test.tsv.

Usage (once per train.tsv/test.tsv):
    python prepare_seed_split_for_models.py \
        --input data/Splits-Rhea/Scaffold/seed_splits/seed0/train.tsv \
        --output data/Splits-Rhea/Scaffold/seed_splits/seed0/prepared/train.tsv

    python prepare_seed_split_for_models.py \
        --input data/Splits-Rhea/Scaffold/seed_splits/seed0/test.tsv \
        --output data/Splits-Rhea/Scaffold/seed_splits/seed0/prepared/test.tsv \
        --queries_output data/Splits-Rhea/Scaffold/seed_splits/seed0/prepared/queries.txt
"""

import argparse
import os

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Rename a Rhea seed split to the lowercase MetaNetX-style schema.")
    parser.add_argument("--input", required=True,
                         help="Rhea seed's train.tsv or test.tsv (columns REACTION_ID, REACTION_SMILES, EC_NUMBER)")
    parser.add_argument("--output", required=True, help="Output path for the renamed tsv")
    parser.add_argument("--queries_output", default=None,
                         help="If given, also write reaction_smiles (one per line, same row order) here "
                              "-- only meaningful when --input is a test.tsv")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    required = {"REACTION_ID", "REACTION_SMILES", "EC_NUMBER"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"'{args.input}' is missing expected column(s): {sorted(missing)} "
                          f"(found: {list(df.columns)})")

    df = df.rename(columns={
        "REACTION_ID": "reaction_id",
        "REACTION_SMILES": "reaction_smiles",
        "EC_NUMBER": "ec",
    })
    # Some consumers (e.g. theia's encode_split_data.py) expect a 'rxn' column
    # instead of 'reaction_smiles' -- write both, same convention as
    # data/Splits-DBs/MetaNetX/scripts/generate_seed_split.py.
    df.insert(2, "rxn", df["reaction_smiles"])
    df = df[["reaction_id", "reaction_smiles", "rxn", "ec"]]

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote '{args.output}' ({len(df)} rows)")

    if args.queries_output:
        os.makedirs(os.path.dirname(os.path.abspath(args.queries_output)), exist_ok=True)
        with open(args.queries_output, "w") as f:
            for s in df["reaction_smiles"]:
                f.write(str(s) + "\n")
        print(f"Wrote '{args.queries_output}' ({len(df)} lines, row-order aligned with '{args.output}')")


if __name__ == "__main__":
    main()
