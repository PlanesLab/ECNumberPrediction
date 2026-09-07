#!/usr/bin/env python3
"""
Merge rxnfp_create.py's RXNFP batches with DRFP fingerprints into the
combined 512-dim embedding DB CLAIRE's contrastive trainer consumes.

Rewritten with real CLI args (was hardcoded to a nonexistent
/scratch/jarcagniriv/Case2/... tree). Reads the same
train_reactions_expanded.csv that rxnfp_create.py wrote, so row order/count
stays aligned with the RXNFP batches.
"""

import argparse
import os
import pickle
from glob import glob
from collections import defaultdict

import numpy as np
from natsort import natsorted
from drfp import DrfpEncoder
import pandas as pd


def truncate_ec(ec_string: str) -> str:
    parts = ec_string.split('.')
    return '.'.join(parts[:3]) if len(parts) >= 3 else ec_string


def main() -> None:
    parser = argparse.ArgumentParser(description="Build CLAIRE's combined RXNFP+DRFP embedding DB.")
    parser.add_argument("--db_dir", required=True, help="Directory containing rxnfp_embs_*.pkl and train_reactions_expanded.csv (rxnfp_create.py's --output_dir)")
    parser.add_argument("--smiles_col", default="reaction_smiles")
    parser.add_argument("--ec_col", default="ec")
    args = parser.parse_args()

    base = args.db_dir

    batch_files = natsorted(glob(os.path.join(base, 'rxnfp_embs_*.pkl')))
    print(f"Found {len(batch_files)} RXNFP batch files")

    rxnfp_list = []
    for fn in batch_files:
        print(f" -> Loading {os.path.basename(fn)}")
        with open(fn, 'rb') as f:
            data = pickle.load(f)
        arr = np.asarray(data['rxnfp'])
        if arr.ndim == 2 and arr.shape[1] != 256:
            raise ValueError(f"Unexpected RXNFP dim {arr.shape[1]}, expected 256")
        rxnfp_list.append(arr)
    rxnfp = np.vstack(rxnfp_list)
    N, rxnfp_dim = rxnfp.shape
    if rxnfp_dim != 256:
        raise ValueError(f"RXNFP should have dim 256 but got {rxnfp_dim}")
    print(f"RXNFP loaded: {N} samples x {rxnfp_dim} dims")

    csv_path = os.path.join(base, "train_reactions_expanded.csv")
    df = pd.read_csv(csv_path, dtype=str)
    if len(df) != N:
        raise ValueError(f"CSV rows ({len(df)}) != RXNFP samples ({N})")
    smiles = df[args.smiles_col].astype(str).tolist()
    raw_ecs = df[args.ec_col].astype(str).tolist()
    ec3_labels = [truncate_ec(ec) for ec in raw_ecs]

    print("Generating 256-bit DRFP fingerprints...")
    drfp_list = DrfpEncoder.encode(smiles, n_folded_length=256)
    drfps = np.vstack(drfp_list)
    if drfps.shape != (N, 256):
        raise ValueError(f"DRFP shape {drfps.shape}, expected ({N},256)")
    print(f"DRFP generated: {drfps.shape[0]} samples x {drfps.shape[1]} dims")

    print("Concatenating to 512-dimension embeddings...")
    combined = np.hstack((rxnfp, drfps))
    if combined.shape != (N, 512):
        raise ValueError(f"Combined shape {combined.shape}, expected ({N},512)")

    emb_dict = defaultdict(list)
    for vec, ec3 in zip(combined, ec3_labels):
        emb_dict[ec3].append(vec)
    emb_dict = {ec: np.vstack(arrs) for ec, arrs in emb_dict.items()}

    out_emb = os.path.join(base, 'esm_emb_dict_ec3.pkl')
    with open(out_emb, 'wb') as f:
        pickle.dump(emb_dict, f)
    print(f"Saved embedding dict: {out_emb}")

    out_labels = os.path.join(base, 'labels_train_ec3.pkl')
    with open(out_labels, 'wb') as f:
        pickle.dump(ec3_labels, f)
    print(f"Saved labels list: {out_labels}")

    out_lookup = os.path.join(base, 'lookup_array_ec3.pkl')
    with open(out_lookup, 'wb') as f:
        pickle.dump(combined, f)
    print(f"Saved lookup array: {out_lookup}")


if __name__ == '__main__':
    main()
