#!/bin/bash
###################################################################################
### Run SIMMER
###################################################################################

#SBATCH --job-name=simmer
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jarcagniriv@unav.es
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH -o results/logs/simmer_%j.out

set -euo pipefail

echo "==================================="
echo "===     Load the Packages       ==="
echo "==================================="
echo "$(date)"
module load Anaconda3

conda activate SIMMERenv

# --- Case 1: KEGG ---
python methods/SIMMER/SIMMER_scripts/simmer_input.py \
    --input-path data/KEGG/kegg_reactions_with_split_can.csv \
    --output-path methods/SIMMER/SIMMER_scripts/input/kegg_input_simmer.csv \
    --reaction-id-col "Reaction ID" \
    --reaction-sep "<=>" \
    --sp_col Equation \
    --smiles-col Isomeric_SMILES \
    --include-ec False

python3 methods/SIMMER/SIMMER_code/SIMMER/SIMMER.py \
    -i methods/SIMMER/SIMMER_code/SIMMER/SIMMER_files \
    -o methods/SIMMER/SIMMER_scripts/output/ResultsKEGG \
    -q methods/SIMMER/SIMMER_scripts/input/kegg_input_simmer.csv

python methods/SIMMER/SIMMER_scripts/ec_predictions.py \
    --input_dir methods/SIMMER/SIMMER_scripts/output/ResultsKEGG \
    --output_file results/Case1/KEGG-1.8K/SIMMER.csv

# --- Case 2: MetaNetX ---
python methods/SIMMER/SIMMER_scripts/simmer_input.py \
    --input-path data/Splits-DBs/MetaNetX/test.tsv \
    --output-path methods/SIMMER/SIMMER_scripts/input/metanetx_input_simmer.csv \
    --reaction-id-col reaction_id \
    --sp_col substrates_products \
    --smiles-col reaction_smiles \
    --include-ec False

python methods/SIMMER/SIMMER_scripts/simmer_input.py \
    --input-path data/Splits-DBs/MetaNetX/train.tsv \
    --output-path methods/SIMMER/SIMMER_code/SIMMER/SIMMER_files_metanetx/chem_data/metanetx_reactions.csv \
    --reaction-id-col reaction_id \
    --sp_col substrates_products \
    --smiles-col reaction_smiles \
    --include-ec True \
    --ec-col ec

python methods/SIMMER/SIMMER_scripts/create_SIMMER_db.py
python methods/SIMMER/SIMMER_scripts/ec_permutations.py

python3 methods/SIMMER/SIMMER_code/SIMMER/SIMMER2.py \
    -i methods/SIMMER/SIMMER_code/SIMMER/SIMMER_files_metanetx \
    -o methods/SIMMER/SIMMER_scripts/output/ResultsMetaNetX \
    -q methods/SIMMER/SIMMER_scripts/input/metanetx_input_simmer.csv

python methods/SIMMER/SIMMER_scripts/ec_predictions.py \
    --input_dir methods/SIMMER/SIMMER_scripts/output/ResultsMetaNetX \
    --output_file results/Case2/results-DBs/MetaNetX/SIMMER.csv

# --- Case Study: drugs ---
python methods/SIMMER/SIMMER_scripts/simmer_input.py \
    --input-path data/Drugs/drug_smiles_updated.csv \
    --output-path methods/SIMMER/SIMMER_scripts/input/simmer_drugs_input.csv \
    --reaction-id-col drug \
    --sp_col right_comp \
    --smiles-col reaction_smiles \
    --include-ec False

python3 methods/SIMMER/SIMMER_code/SIMMER/SIMMER.py \
    -i methods/SIMMER/SIMMER_code/SIMMER/SIMMER_files \
    -o methods/SIMMER/SIMMER_scripts/output/ResultsDrugs \
    -q methods/SIMMER/SIMMER_scripts/input/simmer_drugs_input.csv

python methods/SIMMER/SIMMER_scripts/ec_predictions.py \
    --input_dir methods/SIMMER/SIMMER_scripts/output/ResultsDrugs \
    --output_file results/CaseStudy/results/SIMMER.csv
