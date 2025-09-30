#!/usr/bin/env python3
import os
import pickle
import pandas as pd
from rxnfp.transformer_fingerprints import get_default_model_and_tokenizer, RXNBERTFingerprintGenerator

def main():
    # Paths
    input_csv = '/scratch/jarcagniriv/Case2/MetaNetX/train_reactions_expanded.csv'
    rxnfp_dir = '/scratch/jarcagniriv/Case2/CLAIRE/dev/data/pred_metanetx/'
    rxnfp_prefix = 'rxnfp_embs'
    os.makedirs(rxnfp_dir, exist_ok=True)

    # Read and extract
    df = pd.read_csv(input_csv)
    df['ec3'] = df['ec'].astype(str).str.split('.').str[:3].str.join('.')
    rxns = df['rxn'].tolist()

    # Model setup
    model, tokenizer = get_default_model_and_tokenizer()
    rxnfp_gen = RXNBERTFingerprintGenerator(model, tokenizer)
    print(f"Generating RXNFP embeddings for {len(rxns)} reactions in batches...")

    # Batch parameters
    batch_size = 2000
    all_meta = []

    for i in range(0, len(rxns), batch_size):
        batch_rxns = rxns[i : i + batch_size]
        fps = rxnfp_gen.convert_batch(batch_rxns)  # should be (batch_size, 256)
        batch_meta = {
            'start': i,
            'end': min(i + batch_size, len(rxns)) - 1,
            'rxns': batch_rxns,
            'ec3': df['ec3'].tolist()[i : i + batch_size],
            'rxnfp': fps
        }
        out_pkl = os.path.join(
            rxnfp_dir,
            f"{rxnfp_prefix}_{i:05d}_{i+len(batch_rxns)-1:05d}.pkl"
        )
        with open(out_pkl, 'wb') as f:
            pickle.dump(batch_meta, f)
        print(f"Saved batch {i}-{i+len(batch_rxns)-1} to {out_pkl}")

    print("All batches saved.")

if __name__ == '__main__':
    main()

