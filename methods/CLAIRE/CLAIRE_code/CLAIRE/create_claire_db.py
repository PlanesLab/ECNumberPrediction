#!/usr/bin/env python3
import os
import pickle
import numpy as np
from glob import glob
from natsort import natsorted
from drfp import DrfpEncoder
import pandas as pd
from collections import defaultdict

def truncate_ec(ec_string):
    parts = ec_string.split('.')
    return '.'.join(parts[:3]) if len(parts) >= 3 else ec_string

def main():
    base = '/scratch/jarcagniriv/Case2/CLAIRE/dev/data/pred_metanetx'

    # 1) Merge RXNFP batches (should produce N x 256)
    batch_pattern = os.path.join(base, 'rxnfp_embs_*.pkl')
    batch_files = natsorted(glob(batch_pattern))
    print(f"Found {len(batch_files)} RXNFP batch files")

    rxnfp_list = []
    for fn in batch_files:
        print(f" → Loading {os.path.basename(fn)}")
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
    print(f"RXNFP loaded: {N} samples × {rxnfp_dim} dims")

    # 2) Load reactions and labels
    csv_path = '/scratch/jarcagniriv/Case2/MetaNetX/train_reactions_expanded.csv'
    df = pd.read_csv(csv_path)
    if len(df) != N:
        raise ValueError(f"CSV rows ({len(df)}) != RXNFP samples ({N})")
    smiles = df['rxn'].astype(str).tolist()
    raw_ecs = df['ec'].astype(str).tolist()
    ec3_labels = [truncate_ec(ec) for ec in raw_ecs]

    # 3) Generate 256-bit DRFPs via static encode method
    print("Generating 256-bit DRFP fingerprints…")
    # DrfpEncoder.encode accepts list or single SMILES, returns list of numpy arrays
    drfp_list = DrfpEncoder.encode(smiles, n_folded_length=256)
    drfps = np.vstack(drfp_list)
    if drfps.shape != (N, 256):
        raise ValueError(f"DRFP shape {drfps.shape}, expected ({N},256)")
    print(f"DRFP generated: {drfps.shape[0]} samples × {drfps.shape[1]} dims")

    # 4) Concatenate to combined 512-dimension fingerprint
    print("Concatenating to 512-dimension embeddings…")
    combined = np.hstack((rxnfp, drfps))
    if combined.shape != (N, 512):
        raise ValueError(f"Combined shape {combined.shape}, expected ({N},512)")

    # 5) Group by EC3
    emb_dict = defaultdict(list)
    for vec, ec3 in zip(combined, ec3_labels):
        emb_dict[ec3].append(vec)
    emb_dict = {ec: np.vstack(arrs) for ec, arrs in emb_dict.items()}

    # 6) Save outputs
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
