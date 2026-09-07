"""
Build a SIMMER2-compatible reference database (chem_data/ + prot_data/) from
a no-cofactor Rhea Stratified training split.

Rhea's REACTION_ID is not unique per row (the same base reaction can carry
several different EC labels, one per row -- e.g. broad-specificity kinases).
SIMMER's dictionaries are keyed by the "reaction" id, so a bare REACTION_ID
would silently collide and drop rows. We therefore key each DB entry by
"{REACTION_ID}_{EC_NUMBER}", mirroring how the upstream SIMMER data explodes
multi-EC MetaCyc rows.

Usage:
    python build_rhea_nocofactor_db.py \
        --input ../../../data/Splits-Rhea/Stratified/train_nocofactor.tsv \
        --output-dir .../SIMMER_files_rhea_stratified_nocofactor
"""

import argparse
import os
import pickle

import numpy as np
import pandas as pd
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs


def create_fingerprints(df):
    fps = []
    for _, row in df.iterrows():
        left, right = row["REACTION_SMILES"].split(">>")
        rxn = rdChemReactions.ReactionFromSmarts(f"{left.strip()}>>{right.strip()}")
        fp = rdChemReactions.CreateDifferenceFingerprintForReaction(rxn)
        fps.append(fp)
    return fps


def create_tanimoto_matrix(fps):
    n = len(fps)
    matrix = np.zeros((n, n), dtype=np.float64)
    for i, fp in enumerate(fps):
        matrix[i, :] = DataStructs.BulkTanimotoSimilarity(fp, fps)
        if i % 2000 == 0:
            print(f"Processed row {i} of {n}", flush=True)
    return matrix


def take_ES_walk(ec_cat, ec_df, level):
    tally = 0
    running_tally = []
    for ec in ec_df[level]:
        if ec == ec_cat:
            tally += 1
        elif ec in ["NIL", "DM"]:
            pass
        else:
            tally -= 1
        running_tally.append(tally)
    return running_tally


def compute_perm_score(ec_df, level, ec_cat):
    permuted = ec_df[level].sample(frac=1, replace=False).tolist()
    temp_df = pd.DataFrame({level: permuted})
    es_tally = take_ES_walk(ec_cat, temp_df, level)
    denom = ec_df[level].value_counts().get(ec_cat, 1)
    return max(es_tally) / denom


def build_ec_permutations(ec_series, num_permutations, output_csv):
    ec_df = pd.DataFrame({"EC4": ec_series})
    for lvl in (1, 2, 3):
        ec_df[f"EC{lvl}"] = ec_df["EC4"].apply(lambda x: ".".join(x.split(".")[:lvl]))
    reps = {}
    for lvl in (1, 2, 3, 4):
        col = f"EC{lvl}"
        valid = ec_df[~ec_df[col].isin(["NIL", "DM"])]
        reps[col] = valid[col].mode().iloc[0] if not valid.empty else None

    rows = []
    for _ in range(num_permutations):
        rows.append({col: compute_perm_score(ec_df, col, rep) if rep is not None else None
                     for col, rep in reps.items()})
    perm_df = pd.DataFrame(rows)
    perm_df.to_csv(output_csv, index=False)
    print(f"Permutation scores saved to {output_csv}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True, help="SIMMER_files_* directory")
    parser.add_argument("--num-permutations", type=int, default=1000)
    args = parser.parse_args()

    chem_dir = os.path.join(args.output_dir, "chem_data")
    prot_dir = os.path.join(args.output_dir, "prot_data")
    os.makedirs(chem_dir, exist_ok=True)
    os.makedirs(prot_dir, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    df["EC_NUMBER"] = df["EC_NUMBER"].astype(str).str.strip()
    df = df[df["EC_NUMBER"].apply(lambda x: x and x.lower() != "nan")].reset_index(drop=True)
    # Hyphen (not underscore) keeps the id intact through ec_predictions.py's
    # `filename.split('_')[0]` reaction-name extraction downstream.
    df["reaction"] = df["REACTION_ID"].astype(str) + "-" + df["EC_NUMBER"]
    assert df["reaction"].is_unique, "DB reaction ids are not unique after REACTION_ID+EC composition"

    print(f"Building DB from {len(df)} reactions...")
    fps = create_fingerprints(df)
    with open(os.path.join(chem_dir, "MC_rxn_fps.p"), "wb") as f:
        pickle.dump(fps, f)

    rxn_to_ec = dict(zip(df["reaction"], df["EC_NUMBER"]))
    with open(os.path.join(chem_dir, "MC_rxn_ec_dict.p"), "wb") as f:
        pickle.dump(rxn_to_ec, f)

    tanimoto = create_tanimoto_matrix(fps)
    pd.DataFrame(tanimoto).to_csv(os.path.join(chem_dir, "MC_rxn_tanimoto_matrix.csv"), index=False)

    reaction_to_index = {rxn: idx for idx, rxn in enumerate(df["reaction"])}
    with open(os.path.join(chem_dir, "MC_rxn_to_matrix_index_dict.p"), "wb") as f:
        pickle.dump(reaction_to_index, f)

    # Hardcoded filename expected by SIMMER2.py's main(); content is loaded but
    # unused downstream beyond existing, so we just keep the reaction table.
    df.to_csv(os.path.join(chem_dir, "metanetx_reactions.csv"), index=False)

    build_ec_permutations(df["EC_NUMBER"], args.num_permutations,
                           os.path.join(chem_dir, "ec_perm.csv"))

    with open(os.path.join(prot_dir, "prot_dict.p"), "wb") as f:
        pickle.dump({}, f)

    print("Done.")
    print(f"DB size (for SIMMER2.py reference-size check): {len(df)}")


if __name__ == "__main__":
    main()
