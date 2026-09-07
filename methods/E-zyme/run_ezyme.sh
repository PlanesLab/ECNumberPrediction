#!/bin/bash
# -----------------------------------------------------------------------------
# Author: jarcagniriv
# Script: run_ezyme.sh
# Description: E-zyme 1 and 2 web scraping and result extraction for EC number prediction.
# -----------------------------------------------------------------------------
set -euo pipefail

# --- Case 1: KEGG substrate-product pairs ---
# SEED selects which seed's resampled KEGG-1.8K test set to use (see run_theia.sh for the same
# convention); defaults to the original fixed test set for backward compatibility. Scrapes go
# into a SHARED cache dir across all seeds (seed pools are resampled from the same ~7.8K KEGG
# pool and overlap heavily) -- ezyme_webscrapping.py already skips pairs it's already scraped,
# and --reaction_ids_filter restricts each seed's extracted CSV to just its own subset.
SEED="${SEED:-42}"
if [ -d "data/Subsets/KEGG/seed_pools/seed${SEED}" ]; then
    KEGG_TEST_CSV="data/Subsets/KEGG/seed_pools/seed${SEED}/kegg_reactions_current_test.csv"
    KEGG_OUT_DIR="results/Case1/seed_runs/seed${SEED}/KEGG-1.8K"
else
    KEGG_TEST_CSV="data/Subsets/KEGG/kegg_reactions_current_test.csv"
    KEGG_OUT_DIR="results/Case1/KEGG-1.8K"
fi
mkdir -p "$KEGG_OUT_DIR"

python methods/E-zyme/ezyme_scripts/derive_sp_pairs.py \
    --input "$KEGG_TEST_CSV" \
    --output "methods/E-zyme/output/kegg_sp_pairs_seed${SEED}.csv"

python methods/E-zyme/ezyme_scripts/ezyme_webscrapping.py \
    -i "methods/E-zyme/output/kegg_sp_pairs_seed${SEED}.csv" \
    -o methods/E-zyme/output/outputKEGG_pool \
    --reaction_id_col "Reaction ID" \
    --reactant_col "Reactants" \
    --product_col "Products" \
    --delimiter ","

# After E-zyme server returns results, extract EC numbers (restricted to this seed's subset):
python methods/E-zyme/ezyme_scripts/get_ezyme_results.py \
    --input_dir methods/E-zyme/output/outputKEGG_pool \
    --reaction_ids_filter "$KEGG_TEST_CSV" \
    --output_file "$KEGG_OUT_DIR/E-zyme.csv"

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
