#!/bin/bash
###################################################################################
# Author: Josefina Arcagni
# Date: 9/9/2025
# Description: Run SelenzymeRF for Case 1, Case 2, and Case Study.
###################################################################################
set -euo pipefail

# --- Case 1: KEGG ---
echo "Running SelenzymeRF for Case 1"
bash methods/SelenzymeRF/SelenzymeRF_code/start_server.sh &
SERVER_PID=$!

# Wait for server to start
sleep 300

python3 methods/SelenzymeRF/SelenzymeRF_scripts/query_selenzyme.py \
    --csv_file data/KEGG/kegg_reactions_current_test.csv \
    --results_folder results/SelenzymeRF/case1 \
    --server_url http://localhost:5000 \
    --reaction_name_column "Reaction ID" \
    --reaction_smiles_column Isomeric_SMILES

python3 methods/SelenzymeRF/SelenzymeRF_scripts/get_results.py \
    --folder results/SelenzymeRF/case1 \
    --ec_column "EC Number" \
    --output results/Case1/KEGG-1.8K/SelenzymeRF.csv

kill "$SERVER_PID" 2>/dev/null || true

# --- Case Study ---
echo "Running SelenzymeRF for Case Study"
bash methods/SelenzymeRF/SelenzymeRF_code/start_server.sh &
SERVER_PID=$!

sleep 300

python3 methods/SelenzymeRF/SelenzymeRF_scripts/query_selenzyme.py \
    --csv_file data/Drugs/drug_smiles_updated.csv \
    --results_folder results/SelenzymeRF/casestudy \
    --server_url http://localhost:5000 \
    --reaction_name_column drug \
    --reaction_smiles_column reaction_smiles

python3 methods/SelenzymeRF/SelenzymeRF_scripts/get_results.py \
    --folder results/SelenzymeRF/casestudy \
    --ec_column "EC Number" \
    --output results/CaseStudy/results/SelenzymeRF.csv

kill "$SERVER_PID" 2>/dev/null || true

# --- Case 2: replace SelenzymeRF database with custom training set ---
echo "Running SelenzymeRF for Case 2"
python3 methods/SelenzymeRF/SelenzymeRF_scripts/generate_selenzyme_db.py

bash methods/SelenzymeRF/SelenzymeRF_code/start_server.sh &
SERVER_PID=$!

sleep 300

python3 methods/SelenzymeRF/SelenzymeRF_scripts/query_selenzyme.py \
    --csv_file data/Rhea/test_reactions.tsv \
    --results_folder results/SelenzymeRF/case2 \
    --server_url http://localhost:5000 \
    --reaction_name_column reaction_id \
    --reaction_smiles_column reaction_smiles

python3 methods/SelenzymeRF/SelenzymeRF_scripts/get_results.py \
    --folder results/SelenzymeRF/case2 \
    --ec_column "EC Number" \
    --output results/Case2/results-DBs/Rhea/SelenzymeRF.csv

kill "$SERVER_PID" 2>/dev/null || true
