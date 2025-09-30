"""
Script: generate_selenzyme_db.py
Author: Josefina Arcagni
Date: 2025-09-09
Description: This script generates a new Selenzyme database by filtering out reactions present in the test set.
"""

import pandas as pd
import os

# File paths
test_reactions_file = "/data/MetaNetX/test_reactions.tsv"

input_files = {
    "reac_prop.tsv": "SelenzymeRF/Db_old/reac_prop.tsv",
    "reac_seqs.tsv": "SelenzymeRF/Db_old/reac_seqs.tsv",
    "reac_smi.csv":  "/SelenzymeRF/Db_old/reac_smi.csv"
}

# Output directory
output_dir = "/scratch/jarcagniriv/Case2/SelenzymeRF/DB_new"
os.makedirs(output_dir, exist_ok=True)

# Read the test reactions file (assumed to be TSV)
test_df = pd.read_csv(test_reactions_file, sep='\t')
# Create a set of reaction IDs to filter out
reaction_ids_to_remove = set(test_df["reaction_id"])

# Process each input file
for filename, filepath in input_files.items():
    # Determine the separator based on the file extension
    sep = "\t" if filename.endswith(".tsv") else ","
    
    # Read the file
    df = pd.read_csv(filepath, sep=sep)
    
    # Get the first column name (assuming the reaction code is in the first column)
    first_col = df.columns[0]
    
    # Filter out rows where the value in the first column is in the reaction_ids_to_remove set
    filtered_df = df[~df[first_col].isin(reaction_ids_to_remove)]
    
    # Define the output path
    out_path = os.path.join(output_dir, filename)
    
    # Write the filtered DataFrame to file with the same separator as the input
    filtered_df.to_csv(out_path, sep=sep, index=False)
    
    print(f"Processed {filename} -> {out_path}")

print("All files processed successfully.")
