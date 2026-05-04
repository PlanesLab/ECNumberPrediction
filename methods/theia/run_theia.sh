#!/bin/bash
###################################################################################
# Author: Josefina Arcagni
# Date: 9/9/2025
# Description: Run Theia for Case 1, Case 2, and Case Study.
###################################################################################

#SBATCH --job-name=theia_GPU
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jarcagniriv@unav.es
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH -o results/logs/theia_%j.out

set -euo pipefail

module load Anaconda3
conda activate theia

python -c 'import sys; print(sys.version_info[:])'
echo "$PATH"

# --- Case 1: KEGG ---
python methods/theia/theia_scripts/query_theia.py \
    --query_file "data/Subsets/KEGG/8:2KEGGTest_canonicalized.txt" \
    --reaction_ids_file data/Subsets/KEGG/kegg_reactions_current_test.csv \
    --reaction_id_column "Reaction ID" \
    --output_file results/Case1/KEGG-1.8K/theia.csv

# --- Case 2: MetaNetX ---
# Train a custom Theia model on the MetaNetX training split
python methods/theia/theia_code/theia/generateDB/encode_split_data.py \
    --path data/Splits-DBs/MetaNetX/train.tsv \
    --output-path methods/theia/theia_code/theia/DB

bash methods/theia/theia_code/theia/train_all.sh

python methods/theia/theia_scripts/query_theia.py \
    --query_file data/Splits-DBs/MetaNetX/queries.txt \
    --reaction_ids_file data/Splits-DBs/MetaNetX/test.tsv \
    --reaction_id_column reaction_id \
    --model DB.123 \
    --output_file results/Case2/results-DBs/MetaNetX/theia.csv

# --- Case Study ---
python methods/theia/theia_scripts/query_theia.py \
    --query_file data/Drugs/reaction_smiles_can.txt \
    --reaction_ids_file data/Drugs/drug_smiles_updated.csv \
    --reaction_id_column drug \
    --output_file results/CaseStudy/results/theia.csv
