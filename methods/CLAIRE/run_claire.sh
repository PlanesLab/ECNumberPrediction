#!/bin/bash
###################################################################################
# Author: Josefina Arcagni
# Date: 9/9/2025
# Description: Run CLAIRE for Case 1, Case 2, and Case Study.
###################################################################################

#SBATCH --job-name=claire_GPU
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jarcagniriv@unav.es
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH -o results/logs/claire_%j.out

set -euo pipefail

module load Anaconda3
python -c 'import sys; print(sys.version_info[:])'
echo "$PATH"

CLAIRE_CODE="methods/CLAIRE/CLAIRE_code/CLAIRE"
CLAIRE_SCRIPTS="methods/CLAIRE/CLAIRE_scripts"

# --- Case 1: KEGG ---
echo "==================================="
echo "===  Generate DRFP Fingerprints  ==="
echo "==================================="
echo "$(date)"

QUERY_PATH="data/Subsets/KEGG/8:2KEGGTest_canonicalized.txt"
FPS_PATH="$CLAIRE_CODE/fps/my_rxn_fps.pkl"
TEST_DATA_OUT="$CLAIRE_CODE/Results/claire_test_kegg.npz"

conda activate claire
drfp "$QUERY_PATH" "$FPS_PATH" -d 256
conda deactivate

conda activate rxnfp-env
python "$CLAIRE_SCRIPTS/create_fps.py" \
    --query_path "$QUERY_PATH" \
    --fps_path "$FPS_PATH" \
    --test_out "$TEST_DATA_OUT"
conda deactivate

conda activate claire
python "$CLAIRE_SCRIPTS/query_claire.py" \
    --test_data_path "$TEST_DATA_OUT" \
    --train_data_path "$CLAIRE_CODE/dev/data/model_lookup_train.pkl" \
    --train_labels_path "$CLAIRE_CODE/dev/data/pred_rxn_EC123/labels_train_ec3.pkl" \
    --test_csv_path data/Subsets/KEGG/kegg_reactions_with_split_can.csv \
    --reaction_id_col "Reaction ID" \
    --model_path "$CLAIRE_CODE/dev/results/model/pred_rxn_EC123/layer5_node1280_triplet2000_final.pth" \
    --gmm_path "$CLAIRE_CODE/dev/gmm/gmm_ensumble.pkl"

python "$CLAIRE_SCRIPTS/get_ec_predictions.py" \
    --input_path "$CLAIRE_CODE/Results/results_kegg20.csv" \
    --output_path results/Case1/KEGG-1.8K/CLAIRE.csv

# --- Case 2: MetaNetX ---
echo "==================================="
echo "===  Train CLAIRE model  ==="
echo "==================================="
echo "$(date)"

conda activate claire
python "$CLAIRE_CODE/dev/training/train-pred_rxn_EC.py"

QUERY_PATH="data/Splits-DBs/MetaNetX/queries.txt"
FPS_PATH="$CLAIRE_CODE/fps/my_rxn_fps.pkl"
TEST_DATA_OUT="$CLAIRE_CODE/Results/claire_test_metanetx.npz"

drfp "$QUERY_PATH" "$FPS_PATH" -d 256
conda deactivate

conda activate rxnfp-env
python "$CLAIRE_SCRIPTS/create_fps.py" \
    --query_path "$QUERY_PATH" \
    --fps_path "$FPS_PATH" \
    --test_out "$TEST_DATA_OUT"
conda deactivate

conda activate claire
python "$CLAIRE_SCRIPTS/query_claire.py" \
    --test_data_path "$TEST_DATA_OUT" \
    --train_data_path "$CLAIRE_CODE/dev/data/pred_metanetx/lookup_array_ec3.pkl" \
    --train_labels_path "$CLAIRE_CODE/dev/data/pred_metanetx/labels_train_ec3.pkl" \
    --test_csv_path data/Splits-DBs/MetaNetX/test.tsv \
    --reaction_id_col reaction_id \
    --model_path "$CLAIRE_CODE/dev/results/model/pred_rxn_metanetx/train_final.pth" \
    --gmm_path "$CLAIRE_CODE/dev/gmm/gmm_ensumble.pkl"

python "$CLAIRE_SCRIPTS/get_ec_predictions.py" \
    --input_path "$CLAIRE_CODE/Results/results_metanetx.csv" \
    --output_path results/Case2/results-DBs/MetaNetX/CLAIRE.csv

# --- Case Study ---
QUERY_PATH="data/Drugs/reaction_smiles_can.txt"
FPS_PATH="$CLAIRE_CODE/fps/my_rxn_fps.pkl"
TEST_DATA_OUT="$CLAIRE_CODE/Results/claire_test_drugs.npz"

drfp "$QUERY_PATH" "$FPS_PATH" -d 256
conda deactivate

conda activate rxnfp-env
python "$CLAIRE_SCRIPTS/create_fps.py" \
    --query_path "$QUERY_PATH" \
    --fps_path "$FPS_PATH" \
    --test_out "$TEST_DATA_OUT"
conda deactivate

conda activate claire
python "$CLAIRE_SCRIPTS/query_claire.py" \
    --test_data_path "$TEST_DATA_OUT" \
    --train_data_path "$CLAIRE_CODE/dev/data/model_lookup_train.pkl" \
    --train_labels_path "$CLAIRE_CODE/dev/data/pred_rxn_EC123/labels_train_ec3.pkl" \
    --test_csv_path data/Drugs/drug_smiles_updated.csv \
    --reaction_id_col drug \
    --model_path "$CLAIRE_CODE/dev/results/model/pred_rxn_EC123/layer5_node1280_triplet2000_final.pth" \
    --gmm_path "$CLAIRE_CODE/dev/gmm/gmm_ensumble.pkl"

python "$CLAIRE_SCRIPTS/get_ec_predictions.py" \
    --input_path "$CLAIRE_CODE/Results/results_drugs.csv" \
    --output_path results/CaseStudy/results/CLAIRE.csv
