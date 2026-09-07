"""
Builds ground-truth EC labels for the Rhea-500 Case 1 test set.

Matches data/Subsets/Rhea/reaction_smiles.txt (500 raw SMILES, no IDs) against
data/Splits-Rhea/master.tsv by canonical SMILES to recover reaction ID and EC.
Matches on both the forward and side-swapped canonical SMILES, since master.tsv
keeps only one direction (LR) per transformation.

Writes data/Subsets/Rhea/rhea500_ground_truth.csv (Reaction ID, Isomeric_SMILES,
EC Number), row-order aligned to reaction_smiles.txt. Unmatched rows get an
empty EC Number rather than being dropped, to keep positional alignment intact.
"""

import argparse
import os
import sys

import pandas as pd
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'Preprocessing'))
from canonicalize_rxn_SMILES import canonicalize_reaction_smiles  # noqa: E402

MASTER_POOL = "data/Splits-Rhea/master.tsv"
RHEA500_SMILES = "data/Subsets/Rhea/reaction_smiles.txt"
OUTPUT_CSV = "data/Subsets/Rhea/rhea500_ground_truth.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description="Recover Rhea-500 ground truth via SMILES match.")
    parser.add_argument("--master_pool", default=MASTER_POOL)
    parser.add_argument("--rhea500_smiles", default=RHEA500_SMILES)
    parser.add_argument("--output", default=OUTPUT_CSV)
    args = parser.parse_args()

    master = pd.read_csv(args.master_pool, sep='\t', dtype=str)
    print(f"Master pool: {len(master)} reactions")
    master['canon'] = master['REACTION_SMILES'].astype(str).apply(canonicalize_reaction_smiles)
    lookup = {}
    for _, row in master.iterrows():
        canon = row['canon']
        if canon is None:
            continue
        lhs, _, rhs = canon.partition(">>")
        value = (row['REACTION_ID'], row['EC_NUMBER'])
        lookup.setdefault(canon, value)
        lookup.setdefault(f"{rhs}>>{lhs}", value)  # direction-invariant: reverse also resolves
    print(f"Unique canonical SMILES in master pool (both directions): {len(lookup)}")

    with open(args.rhea500_smiles) as f:
        rhea_smiles = [line.strip() for line in f if line.strip()]
    print(f"Rhea-500 reactions: {len(rhea_smiles)}")

    records = []
    n_matched = 0
    for i, smi in enumerate(rhea_smiles):
        canon = canonicalize_reaction_smiles(smi)
        rxn_id, ec = lookup.get(canon, (None, None))
        if rxn_id is not None:
            n_matched += 1
        else:
            rxn_id = f"RHEA500_{i}"
            ec = ""
        records.append({"Reaction ID": rxn_id, "Isomeric_SMILES": smi, "EC Number": ec})

    out_df = pd.DataFrame(records)
    out_df.to_csv(args.output, index=False)
    print(f"Matched {n_matched}/{len(rhea_smiles)} reactions against the master pool")
    print(f"Wrote {len(out_df)} rows to '{args.output}' (row-order aligned with reaction_smiles.txt)")


if __name__ == "__main__":
    main()
