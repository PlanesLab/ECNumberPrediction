#!/bin/bash
# -----------------------------------------------------------------------------
# Author: jarcagniriv
# Script: run_BridgIT.sh
# Description: Prepare BridgIT web server input and extract EC numbers from results.
# NOTE: BridgIT requires uploading input ZIPs to https://lcsb-databases.epfl.ch/Bridgit
#       (user account required). Download the result ZIPs before running get_results.py.
# -----------------------------------------------------------------------------
set -euo pipefail

# --- Case 1: KEGG ---
python methods/BridgIT/bridgIT_scripts/bridgit_input.py \
    --input_file data/KEGG/kegg_reactions_current_test.csv \
    --output_file methods/BridgIT/input/reduced_inputs_KEGG/systemfile.txt \
    --molfile_folder methods/BridgIT/input/reduced_inputs_KEGG/molfiles \
    --equation_column Equation

# For large inputs, split into smaller batches first:
#   python methods/BridgIT/bridgIT_scripts/bridgit_reduced_input.py

# Upload the ZIP files to https://lcsb-databases.epfl.ch/Bridgit, then extract:
python methods/BridgIT/bridgIT_scripts/get_results.py \
    --bridgit_dir methods/BridgIT/input/reduced_inputs_KEGG \
    --input_csv data/KEGG/kegg_reactions_current_test.csv \
    --output results/Case1/KEGG-1.8K/BridgIT.csv \
    --csv_reaction_col "Reaction ID" \
    --csv_ec_col "EC Number"

# --- Case Study ---
python methods/BridgIT/bridgIT_scripts/bridgit_input_drugs.py

python methods/BridgIT/bridgIT_scripts/get_results.py \
    --bridgit_dir methods/BridgIT/input/input_casestudy \
    --input_csv data/Drugs/drug_smiles_updated.csv \
    --output results/CaseStudy/results/BridgIT.csv \
    --csv_reaction_col drug \
    --csv_ec_col ec
