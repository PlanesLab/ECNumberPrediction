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
# SEED selects which seed's resampled KEGG-1.8K test set to use (see run_theia.sh for the same
# convention); defaults to the original fixed test set for backward compatibility. Each seed
# needs its own manual upload/download round at https://lcsb-databases.epfl.ch/Bridgit.
SEED="${SEED:-42}"
if [ -d "data/Subsets/KEGG/seed_pools/seed${SEED}" ]; then
    KEGG_TEST_CSV="data/Subsets/KEGG/seed_pools/seed${SEED}/kegg_reactions_current_test.csv"
    KEGG_OUT_DIR="results/Case1/seed_runs/seed${SEED}/KEGG-1.8K"
else
    KEGG_TEST_CSV="data/Subsets/KEGG/kegg_reactions_current_test.csv"
    KEGG_OUT_DIR="results/Case1/KEGG-1.8K"
fi
mkdir -p "$KEGG_OUT_DIR"
BRIDGIT_IN_DIR="methods/BridgIT/input/reduced_inputs_KEGG_seed${SEED}"

python methods/BridgIT/bridgIT_scripts/bridgit_input.py \
    --input_file "$KEGG_TEST_CSV" \
    --output_file "$BRIDGIT_IN_DIR/systemfile.txt" \
    --molfile_folder "$BRIDGIT_IN_DIR/molfiles" \
    --equation_column Equation

# For large inputs, split into smaller batches first:
#   python methods/BridgIT/bridgIT_scripts/bridgit_reduced_input.py

# Upload the ZIP files to https://lcsb-databases.epfl.ch/Bridgit, then extract:
python methods/BridgIT/bridgIT_scripts/get_results.py \
    --bridgit_dir "$BRIDGIT_IN_DIR" \
    --input_csv "$KEGG_TEST_CSV" \
    --output "$KEGG_OUT_DIR/BridgIT.csv" \
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
