# ==============================================================================
# Author: Josefina Arcagni
# Email: jarcagniriv@unav.es
# Date: 2025-02-07
# Script: bridgit_reduced_input.py
# Description: Splits input data into smaller batches, cleans reaction equations by
#              removing variable stoichiometry labels (n, m, etc.) and all spaces.
# ==============================================================================

import os
import shutil
import math
import zipfile
import re

systemfile_path = "/FormattedReactions.txt"
molfiles_path = "/molfiles"
output_base_path = "/methods/BridgIT/input/reduced_inputs_KEGG"

os.makedirs(output_base_path, exist_ok=True)

with open(systemfile_path, "r") as infile:
    lines = infile.readlines()

header = lines[0:4]
reactions = lines[4:]
num_reactions = len(reactions)
split_size = math.ceil(num_reactions / 20)

def clean_equation(equation):
    # Remove variable stoichiometry labels: (n), (m), (n+1), (m-1), (n-m), (m+1), etc.
    equation = re.sub(r"\(n[+-]?\d*\)", "", equation)
    equation = re.sub(r"\(m[+-]?\d*\)", "", equation)
    # Remove all spaces from the equation
    equation = equation.replace(" ", "")
    return equation

for i in range(20):
    start_idx = i * split_size
    end_idx = min((i + 1) * split_size, num_reactions)
    selected_reactions = []
    molecule_ids = set()
    
    for reaction in reactions[start_idx:end_idx]:
        parts = reaction.strip().split(";")
        if len(parts) > 2:
            parts[2] = clean_equation(parts[2])
        cleaned_line = ";".join(parts) + "\n"
        selected_reactions.append(cleaned_line)
        
        for token in parts[2].split("+"):
            token = token.split("<=>")[0]
            token = token.replace("(", "").replace(")", "").strip()
            if token.isalnum():
                molecule_ids.add(token)
    
    split_folder = os.path.join(output_base_path, f"reducedinput{i + 1}")
    os.makedirs(split_folder, exist_ok=True)
    
    split_molfiles_folder = os.path.join(split_folder, "molfiles")
    os.makedirs(split_molfiles_folder, exist_ok=True)
    
    for mol_id in molecule_ids:
        mol_file = f"{mol_id}.mol"
        source_path = os.path.join(molfiles_path, mol_file)
        destination_path = os.path.join(split_molfiles_folder, mol_file)
        if os.path.exists(source_path):
            shutil.copy(source_path, destination_path)
    
    output_file_path = os.path.join(split_folder, "reduced_systemfile.txt")
    with open(output_file_path, "w") as outfile:
        outfile.write("COMPOUNDS\nENTRY\nreactionsS\n")  # Add the requested lines at the beginning
        outfile.writelines(header)
        outfile.writelines(selected_reactions)
    
    zip_file_path = os.path.join(output_base_path, f"reducedinput{i + 1}.zip")
    with zipfile.ZipFile(zip_file_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(split_folder):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, split_folder)
                zipf.write(file_path, arcname)
    
    shutil.rmtree(split_folder)
    print(f"Reduced system file saved to: {zip_file_path}")
