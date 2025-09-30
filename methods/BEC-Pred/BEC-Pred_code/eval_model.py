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
args = parser.parse_args()

logger = logging.getLogger(__name__)

# Path to the trained model directory
train_model_path = args.model_path
model = SmilesClassificationModel("bert", train_model_path, use_cuda=torch.cuda.is_available())

# Load the queries from the text file (one query per row, no header)
with open("/scratch/jarcagniriv/FinalCodes/Case1/KEGG/8:2KEGGTest_canonicalized.txt", "r") as file:
    queries = [line.strip() for line in file if line.strip()]

# Load the reaction IDs from the CSV file from the "Reaction ID" column
reaction_ids_df = pd.read_csv("/scratch/jarcagniriv/FinalCodes/Case1/KEGG/kegg_reactions_current_test.csv")
reaction_ids = reaction_ids_df["Reaction ID"].tolist()

# Predict for the loaded queries
predictions = model.predict(queries)

print("Number of queries:", len(queries))
print("Number of reaction IDs:", len(reaction_ids))

# Create a DataFrame with reaction IDs and predictions
results_df = pd.DataFrame({
    "Reaction ID": reaction_ids,
    "Prediction": predictions[0]
})
print(args.output_csv)
results_df.to_csv(args.output_csv, index=False)
# Print the first few rows of the DataFrame to verify the output
print(results_df.head())
