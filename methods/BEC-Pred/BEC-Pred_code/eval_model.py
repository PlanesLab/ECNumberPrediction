import os
import numpy as np
import pandas as pd
import torch
import logging
import pkg_resources
from rxnfp.models import SmilesClassificationModel
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--model_path", type=str, required=True, help="Path to finetuned model directory")
parser.add_argument("--output_csv", type=str, required=True, help="Path to output CSV file")
parser.add_argument("--queries", type=str, required=True, help="Path to text file with one reaction SMILES per line")
# Optional: without these, results_df carries no reaction ID at all (just a bare Prediction
# column) -- join_results.py needs an ID column to merge on, so any caller that will later join
# these results against ground truth must pass both of these (queries.txt is defined to be
# row-order-aligned with this file, same convention as the reaction_ids_file/queries.txt pairing
# already used by query_theia.py).
parser.add_argument("--reaction_ids_file", type=str, default=None,
                     help="Optional TSV/CSV whose --reaction_id_column, read in the same row order as --queries, is attached to the output as a 'reaction_id' column.")
parser.add_argument("--reaction_id_column", type=str, default="reaction_id")
args = parser.parse_args()

logger = logging.getLogger(__name__)

train_model_path = args.model_path
model = SmilesClassificationModel("bert", train_model_path, use_cuda=torch.cuda.is_available())

with open(args.queries, "r") as file:
    queries = [line.strip() for line in file if line.strip()]

predictions = model.predict(queries)

print("Number of queries:", len(queries))

results_data = {"Prediction": predictions[0]}
if args.reaction_ids_file:
    ids_df = pd.read_csv(args.reaction_ids_file, sep=None, engine="python", dtype=str)
    if len(ids_df) != len(queries):
        raise ValueError(
            f"--reaction_ids_file has {len(ids_df)} rows but --queries has {len(queries)} -- "
            "they must be row-order-aligned (same source split, e.g. queries.txt/test.tsv from prepare_seed_split_for_models.py)."
        )
    results_data = {"reaction_id": ids_df[args.reaction_id_column].tolist(), **results_data}

results_df = pd.DataFrame(results_data)
results_df.to_csv(args.output_csv, index=False)
print(results_df.head())
