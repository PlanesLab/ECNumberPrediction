# ==============================================================================
# Author: Josefina Arcagni
# Email: jarcagniriv@unav.es
# Date: 2025-02-07
# Script: bridgit_input_drugs.py
# Description: Prepares drug reaction data for input into the BridgIT LCSB server.
#              For each reaction, splits the reaction SMILES into individual compounds,
#              assigns unique compound labels (c1, c2, ...), converts each SMILES
#              to a MolFile (using RDKit), and writes a formatted reaction equation.
# ==============================================================================

import os
import csv
import pandas as pd
from rdkit import Chem

# Paths
input_file = '/data/Drugs/drug_smiles_updated.csv'
output_file = '/methods/BridgIT/input/input_casestudy/FormattedReactions.txt'
molfile_folder = '/methods/BridgIT/input/input_casestudy/molfiles'

os.makedirs(molfile_folder, exist_ok=True)

# Global dictionary mapping compound SMILES to invented compound labels (e.g., c1, c2, ...)
compound_map = {}
compound_counter = 1

entries = []

# Read the CSV (assuming comma delimiter)
df = pd.read_csv(input_file)

for index, row in df.iterrows():
    # Use the "drug" column as the reaction identifier (or you can use index)
    reaction_id = row["drug"]
    reaction_smiles = row["reaction_smiles"]

    # Split the reaction SMILES into reactants and products (expecting a single ">>")
    parts = reaction_smiles.split(">>")
    if len(parts) != 2:
        print(f"Reaction {reaction_id} has unexpected format: {reaction_smiles}")
        continue
    left_smiles = parts[0].strip()
    right_smiles = parts[1].strip()

    # Split each side by '.' to get individual compound SMILES
    left_compounds = [s.strip() for s in left_smiles.split('.') if s.strip()]
    right_compounds = [s.strip() for s in right_smiles.split('.') if s.strip()]

    # Process left side compounds
    left_labels = []
    for smi in left_compounds:
        if smi not in compound_map:
            compound_label = f"c{compound_counter}"
            compound_counter += 1
            compound_map[smi] = compound_label
            # Convert SMILES to Mol using RDKit
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                print(f"Error: Could not convert SMILES: {smi} for reaction {reaction_id}")
                continue
            molfile = Chem.MolToMolBlock(mol)
            # Save molfile to disk
            molfile_path = os.path.join(molfile_folder, f"{compound_label}.mol")
            with open(molfile_path, 'w') as f:
                f.write(molfile)
            print(f"Saved molfile for {compound_label} from SMILES: {smi}")
        left_labels.append(compound_map[smi])

    # Process right side compounds
    right_labels = []
    for smi in right_compounds:
        if smi not in compound_map:
            compound_label = f"c{compound_counter}"
            compound_counter += 1
            compound_map[smi] = compound_label
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                print(f"Error: Could not convert SMILES: {smi} for reaction {reaction_id}")
                continue
            molfile = Chem.MolToMolBlock(mol)
            molfile_path = os.path.join(molfile_folder, f"{compound_label}.mol")
            with open(molfile_path, 'w') as f:
                f.write(molfile)
            print(f"Saved molfile for {compound_label} from SMILES: {smi}")
        right_labels.append(compound_map[smi])

    # Create the formatted reaction equation using the compound labels
    formatted_reaction = "+".join(left_labels) + "<=>" + "+".join(right_labels)
    # Format the entry line as: ReactionID;;FormattedReaction;
    entry_line = f"{reaction_id};;{formatted_reaction};"
    entries.append(entry_line)

# Write the formatted reactions to the output file
with open(output_file, 'w') as outfile:
    outfile.write("ENTRY;REACTION;\n")  # header line
    outfile.write("\n".join(entries))

print(f"File saved successfully to {output_file}")
