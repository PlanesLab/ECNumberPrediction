#!/usr/bin/env python3
# ==============================================================================
# Author: Josefina Arcagni (adapted)
# Email: jarcagniriv@unav.es
# Date: 2025-02-07
# Script: generate_becpred_db.py 
# Description:
#   - Processes two KEGG TSV files to generate a database for BEC-Pred.
#   - train_reactions.tsv is processed with all reactions assigned to train.
#   - test_reactions.tsv is processed with all reactions assigned to val.
#   - Both files share the same label mapping (from the first three parts of the EC Number).
#   - Optionally, rows with incomplete EC numbers (if the third dot-separated part is missing or non-numeric)
#     can be removed via a command-line flag.
#   - The label mapping is stored to a CSV file.
# ==============================================================================

import pandas as pd
import random
import argparse

random.seed(41)
label_mapping = {}
current_label = 1

def split_promiscuous_reactions(row):
    """
    If the 'EC_number' field contains multiple EC numbers (separated by pipes '|'),
    return one copy of the row per EC number; otherwise, return the original row in a list.
    """
    if '|' in row['EC_number']:
        ec_numbers = row['EC_number'].split('|')
        new_rows = []
        for ec in ec_numbers:
            new_row = row.copy()
            new_row['EC_number'] = ec.strip()
            new_rows.append(new_row)
        return new_rows
    else:
        return [row]

def split_smiles(reaction_smiles):
    """
    Splits a reaction SMILES string into reactants and products.
    Returns (reactants, products) or (None, None) if splitting fails.
    """
    try:
        reactants, products = reaction_smiles.split('>>')
        return reactants, products
    except ValueError:
        return None, None  

def get_class_label(ec_class_th):
    """
    Returns a numeric label for the EC class (first three parts of an EC Number).
    If the class is new, a new label is assigned.
    """
    global current_label, label_mapping
    if ec_class_th not in label_mapping:
        label_mapping[ec_class_th] = current_label
        current_label += 1
    return label_mapping[ec_class_th]

def process_file(file_path, split_mode, remove_incomplete):
    """
    Processes a TSV file with reaction data.
    
    Parameters:
      file_path (str): Path to the TSV file.
      split_mode (str): If 'random', each valid reaction is assigned a random split
                        (80% train, 20% val). Otherwise, split_mode is used as a fixed
                        split label (e.g., 'val').
      remove_incomplete (bool): If True, skip rows where the EC number is “incomplete” 
                        (i.e., when splitting by '.' the third element is missing or non-numeric).
                        
    Returns:
      A list of dictionaries representing the processed reactions.
    """
    df = pd.read_csv(file_path, sep='\t')
    # Ensure the EC Number column is treated as a string
    df['EC_number'] = df['EC_number'].astype(str)
    # Remove rows with empty EC numbers or those equal to "nan" (case-insensitive)
    df = df[(df['EC_number'].str.strip() != '') & 
            (df['EC_number'].str.strip().str.lower() != 'nan')]
    
    new_rows = []
    
    for _, row in df.iterrows():
        # Skip rows if EC Number is '-' or if the SMILES string is missing
        if row['EC_number'] == '-' or pd.isna(row.get('reaction_smiles', None)):
            continue
        
        # Expand promiscuous reactions (rows with multiple EC numbers separated by pipes)
        for split_row in split_promiscuous_reactions(row):
            ec = split_row['EC_number'].strip()
            # If remove_incomplete is true, skip rows with incomplete EC numbers.
            if remove_incomplete:
                parts = ec.split('.')
                if len(parts) < 3 or not parts[2].strip().isdigit():
                    continue
            
            reactants, products = split_smiles(split_row['reaction_smiles'])
            if reactants is None or products is None or not reactants.strip() or not products.strip():
                continue
            
            # Define the EC class using the first three parts.
            ec_class_th = '.'.join(ec.split('.')[:3])
            class_label = get_class_label(ec_class_th)
            
            if split_mode == 'random':
                split_label = 'train' if random.random() < 0.9 else 'val'
            else:
                split_label = split_mode
            
            new_rows.append({
                'idx': len(new_rows) + 1,
                'rxn': split_row['reaction_smiles'],
                'rxn_class': ec,
                'r': reactants,
                'p': products,
                'rxn_class_th': ec_class_th,
                'class_id': class_label,
                'split': split_label
            })
    
    return new_rows

if __name__ == "__main__":
    # Set up argument parsing.
    parser = argparse.ArgumentParser(
        description="Generate BEC-Pred Database with train/val split from KEGG data."
    )
    parser.add_argument("--remove_incomplete", action="store_true",
                        help="Remove rows with incomplete EC numbers (i.e., missing third part or non-numeric third part)")
    args = parser.parse_args()
    
    file_2018 = '/scratch/jarcagniriv/Case2/MetaNetX/train_reactions.tsv'
    train_val_rows = process_file(file_2018, split_mode='random', remove_incomplete=args.remove_incomplete)
    
    file_newrxns = '/scratch/jarcagniriv/Case2/MetaNetX/test_reactions.tsv'
    newrxns_rows = process_file(file_newrxns, split_mode='val', remove_incomplete=args.remove_incomplete)
    
    df_train_val = pd.DataFrame(train_val_rows)
    df_newrxns = pd.DataFrame(newrxns_rows)
    
    label_data = [{'EC Class': ec_class, 'Assigned Label': label} for ec_class, label in label_mapping.items()]
    df_labels = pd.DataFrame(label_data)
    
    # Save the processed data and label mapping to CSV files.
    output_dir = '/scratch/jarcagniriv/Case2/BEC-Pred/DB/'
    df_train_val.to_csv(output_dir + 'train_metanetx.csv', index=False)
    df_newrxns.to_csv(output_dir + 'test_metanetx.csv', index=False)
    df_labels.to_csv(output_dir + 'ec_class_labels.csv', index=False)
    
    print("Data saved to CSV files successfully.")