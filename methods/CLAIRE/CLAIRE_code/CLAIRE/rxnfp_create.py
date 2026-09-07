#!/usr/bin/env python3
"""
Generate RXNFP embeddings for a MetaNetX train split, batched to pickle
files consumed by create_claire_db.py.

Rewritten with real CLI args (was fully hardcoded to a nonexistent
/scratch/jarcagniriv/Case2/... tree, no --seed). Also expands promiscuous
reactions (EC field with '|'-separated multiple ECs) into one row per EC,
same convention as generate_becpred_db.py's split_promiscuous_reactions,
and writes the expanded frame alongside the embeddings so
create_claire_db.py re-reads the exact same row set/order.
"""

import argparse
import os
import pickle

import pandas as pd
from rxnfp.transformer_fingerprints import get_default_model_and_tokenizer, RXNBERTFingerprintGenerator


def expand_promiscuous(df: pd.DataFrame, smiles_col: str, ec_col: str) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        ec_val = str(row[ec_col])
        for ec in ec_val.split('|'):
            ec = ec.strip()
            if not ec or ec.lower() == 'nan':
                continue
            new_row = row.copy()
            new_row[ec_col] = ec
            rows.append(new_row)
    return pd.DataFrame(rows).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate RXNFP embeddings for CLAIRE Case 2 training.")
    parser.add_argument("--train_file", required=True, help="MetaNetX train.tsv (columns reaction_smiles/rxn, ec)")
    parser.add_argument("--smiles_col", default="reaction_smiles")
    parser.add_argument("--ec_col", default="ec")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--batch_size", type=int, default=2000)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    sep = '\t' if args.train_file.endswith('.tsv') else ','
    df = pd.read_csv(args.train_file, sep=sep, dtype=str)
    df = expand_promiscuous(df, args.smiles_col, args.ec_col)
    df['ec3'] = df[args.ec_col].astype(str).str.split('.').str[:3].str.join('.')

    expanded_csv = os.path.join(args.output_dir, "train_reactions_expanded.csv")
    df.to_csv(expanded_csv, index=False)
    print(f"Expanded {len(df)} rows (promiscuous reactions split), saved to '{expanded_csv}'")

    rxns = df[args.smiles_col].tolist()

    model, tokenizer = get_default_model_and_tokenizer()
    rxnfp_gen = RXNBERTFingerprintGenerator(model, tokenizer)
    print(f"Generating RXNFP embeddings for {len(rxns)} reactions in batches...")

    for i in range(0, len(rxns), args.batch_size):
        batch_rxns = rxns[i: i + args.batch_size]
        fps = rxnfp_gen.convert_batch(batch_rxns)
        batch_meta = {
            'start': i,
            'end': min(i + args.batch_size, len(rxns)) - 1,
            'rxns': batch_rxns,
            'ec3': df['ec3'].tolist()[i: i + args.batch_size],
            'rxnfp': fps,
        }
        out_pkl = os.path.join(args.output_dir, f"rxnfp_embs_{i:05d}_{i + len(batch_rxns) - 1:05d}.pkl")
        with open(out_pkl, 'wb') as f:
            pickle.dump(batch_meta, f)
        print(f"Saved batch {i}-{i + len(batch_rxns) - 1} to {out_pkl}")

    print("All batches saved.")


if __name__ == '__main__':
    main()
