"""
Rebuild master.tsv from raw Rhea distribution files, replacing the previous
SMILES-reversal heuristic (dedup_reversible.py applied after the fact to an
inherited, undocumented master.tsv).

Rhea assigns 4 RHEA IDs to every unique chemical transformation: an
undirected MASTER id, plus LR (forward), RL (reverse) and BI (bidirectional)
variants -- so a naive parse of rhea-reaction-smiles.tsv carries every
transformation twice (once as LR, once as RL). EC numbers are annotated only
once per MASTER id (rhea2ec.tsv, DIRECTION == UN). This script uses
rhea-directions.tsv (Rhea's own master/LR/RL/BI grouping) to collapse each
transformation to exactly one row -- keyed by the MASTER id's EC and the LR
id's SMILES (LR chosen arbitrarily but consistently as "the" direction) --
so no reaction appears twice and no SMILES-matching heuristic is needed.

Masters with more than one distinct EC (e.g. broad-specificity enzymes, 229
of them) are NOT dropped: each is written as multiple rows in master.tsv --
same REACTION_ID and REACTION_SMILES, one row per EC -- so the reaction is
represented once per label rather than once per label removed. These
REACTION_IDs are also listed in master_multi_ec_reactions.tsv so downstream
split scripts can group them together if they want to avoid the same
reaction landing in both train and test under different EC rows.

Masters with no SMILES (LR id has no entry in rhea-reaction-smiles.tsv --
usually generic/macromolecule reactions Rhea encodes without a concrete
SMILES) are excluded and written to master_excluded_no_smiles.tsv instead of
being silently dropped.

Writes data/Splits-Rhea/master.tsv with columns REACTION_ID, REACTION_SMILES,
EC_NUMBER (REACTION_ID is the LR id, matching the existing repo convention --
verified against the previous master.tsv, whose REACTION_ID values are LR ids).
"""

import argparse
import csv
import os
from collections import defaultdict

RAW_DIR = "data/Splits-Rhea/raw"
OUTPUT = "data/Splits-Rhea/master.tsv"


def main() -> None:
    parser = argparse.ArgumentParser(description="Rebuild master.tsv from raw Rhea TSV files.")
    parser.add_argument("--raw_dir", default=RAW_DIR)
    parser.add_argument("--output", default=OUTPUT)
    args = parser.parse_args()

    directions_path = os.path.join(args.raw_dir, "rhea-directions.tsv")
    ec_path = os.path.join(args.raw_dir, "rhea2ec.tsv")
    smiles_path = os.path.join(args.raw_dir, "rhea-reaction-smiles.tsv")

    master_to_lr = {}
    with open(directions_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            master_to_lr[row["RHEA_ID_MASTER"]] = row["RHEA_ID_LR"]
    print(f"Direction groups: {len(master_to_lr)}")

    master_to_ecs = defaultdict(list)
    with open(ec_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            master_to_ecs[row["MASTER_ID"]].append(row["ID"])
    print(f"Masters with an EC annotation: {len(master_to_ecs)}")

    smiles_by_id = {}
    with open(smiles_path) as f:
        for row in csv.reader(f, delimiter="\t"):
            smiles_by_id[row[0]] = row[1]
    print(f"Reaction SMILES rows (LR+RL only): {len(smiles_by_id)}")

    rows, multi_ec_report, no_smiles = [], [], []
    for master_id, ecs in master_to_ecs.items():
        lr_id = master_to_lr.get(master_id)
        unique_ecs = sorted(set(ecs))
        smiles = smiles_by_id.get(lr_id) if lr_id else None
        if smiles is None:
            no_smiles.append({"MASTER_ID": master_id, "LR_ID": lr_id, "EC_NUMBERS": ";".join(unique_ecs)})
            continue
        if len(unique_ecs) > 1:
            multi_ec_report.append({"REACTION_ID": lr_id, "EC_NUMBERS": ";".join(unique_ecs)})
        for ec in unique_ecs:
            rows.append({"REACTION_ID": lr_id, "REACTION_SMILES": smiles, "EC_NUMBER": ec})

    rows.sort(key=lambda r: (int(r["REACTION_ID"]), r["EC_NUMBER"]))

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["REACTION_ID", "REACTION_SMILES", "EC_NUMBER"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} curated rows to '{args.output}' "
          f"(one row per reaction-EC pair; {len(multi_ec_report)} reactions with >1 EC "
          f"contributed more than one row each)")

    excluded_dir = os.path.dirname(args.output)
    multi_ec_path = os.path.join(excluded_dir, "master_multi_ec_reactions.tsv")
    no_smiles_path = os.path.join(excluded_dir, "master_excluded_no_smiles.tsv")

    with open(multi_ec_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["REACTION_ID", "EC_NUMBERS"], delimiter="\t")
        writer.writeheader()
        writer.writerows(multi_ec_report)
    print(f"Wrote {len(multi_ec_report)} multi-EC REACTION_IDs (informational, not excluded) to '{multi_ec_path}'")

    with open(no_smiles_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["MASTER_ID", "LR_ID", "EC_NUMBERS"], delimiter="\t")
        writer.writeheader()
        writer.writerows(no_smiles)
    print(f"Wrote {len(no_smiles)} no-SMILES masters (excluded) to '{no_smiles_path}'")


if __name__ == "__main__":
    main()
