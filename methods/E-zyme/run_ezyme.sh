#!/bin/bash
# -----------------------------------------------------------------------------
# Author: jarcagniriv
# Script: run_ezyme.sh
# Description: E-zyme 1 and 2 web scraping and result extraction for EC number prediction.
# -----------------------------------------------------------------------------
set -euo pipefail

# --- Case 1: KEGG substrate-product pairs ---
python methods/E-zyme/ezyme_scripts/ezyme_webscrapping.py \
    -i data/KEGG/SubstrateProductPairs82.csv \
    -o methods/E-zyme/output/outputKEGG \
    --reaction_id_col "Reaction ID" \
    --reactant_col "Reactants" \
    --product_col "Products" \
    --delimiter ";"

# After E-zyme server returns results, extract EC numbers:
python methods/E-zyme/ezyme_scripts/get_ezyme_results.py \
    --input_dir methods/E-zyme/output/outputKEGG \
    --output_file results/Case1/KEGG-1.8K/E-zyme.csv

# --- Case Study: manually curated drug degradation reactions ---
python methods/E-zyme/ezyme_scripts/ezyme_webscrapping.py \
    -i data/Drugs/sp_pairs_drugs.csv \
    -o methods/E-zyme/output/outputDrugs \
    --reaction_id_col drug \
    --reactant_col Pair1 \
    --product_col Pair2 \
    --delimiter ";"

python methods/E-zyme/ezyme_scripts/get_ezyme_results.py \
    --input_dir methods/E-zyme/output/outputDrugs \
    --output_file results/CaseStudy/results/E-zyme.csv
