#!/bin/bash

###################################################################################
# Author: Josefina Arcagni
# Date: 9/9/2025
# Description: Run CLAIRE for various cases.
###################################################################################


#SBATCH --job-name=claire_GPU_82
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jarcagniriv@unav.es
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --gres=gpu:1      
#SBATCH -o /scratch/jarcagniriv/ECNumberPredictionReview/Claire/logs/bench_%j.out

module load Anaconda3  

python -c 'import sys; print(sys.version_info[:])'

# Optional: Display the PATH for debugging purposes
echo $PATH

################################################### CASE 1 ###########################################################################

echo ===================================
echo ===  Generate DRFP Fingerprints  ===
echo ===================================
echo `date`

QUERY_PATH=/data/KEGG/8:2KEGGTest_canonicalized.txt
FPS_PATH=/CLAIRE_code/CLAIRE/fps/my_rxn_fps.pkl
conda activate claire
drfp "$QUERY_PATH"  "$FPS_PATH" -d 256
conda deactivate 

echo ================================================
echo ======  Generate fps and concatenate  ===========
echo =================================================
echo `date` 

TEST_DATA_OUT= /methods/CLAIRE/CLAIRE_code/CLAIRE/Results/claire_test_kegg.npz
conda activate rxnfp-env
python methods/CLAIRE/CLAIRE_scripts/create_fps.py --query_path "$QUERY_PATH" --fps_path "$FPS_PATH" --test_out "$TEST_DATA_OUT"

echo ===================================
echo ======  Run CLAIRE  ===========
echo ===================================
conda activate claire
python /methods/CLAIRE/CLAIRE_scripts/query_claire.py --test_data_path "$TEST_DATA_OUT" --train_data_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/data/model_lookup_train.pkl 
    / --train_labels_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/data/pred_rxn_EC123/labels_train_ec3.pkl --test_csv_path /data/KEGG/kegg_reactions_with_split_can.csv 
    / --reaction_id_col Reaction ID --model_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/results/model/pred_rxn_EC123/layer5_node1280_triplet2000_final.pth --gmm_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/gmm/gmm_ensumble.pkl

python /methods/CLAIRE/CLAIRE_scripts/get_ec_predictions.py --input_path /methods/CLAIRE/CLAIRE_code/CLAIRE/Results/results_kegg20.csv 
    / --output_path results/Case1/CLAIRE.csv 

################################################### CASE 2 ###########################################################################

echo ===================================
echo ===  Train CLAIRE model  ===
echo ===================================
echo `date`
conda activate claire
python /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/training/train-pred_rxn_EC.py

echo ===================================
echo ===  Generate DRFP Fingerprints  ===
echo ===================================
echo `date`

QUERY_PATH=/data/MetaNetX/queries.txt
FPS_PATH=/CLAIRE_code/CLAIRE/fps/my_rxn_fps.pkl

drfp "$QUERY_PATH"  "$FPS_PATH" -d 256
conda deactivate 

echo ================================================
echo ======  Generate fps and concatenate  ===========
echo =================================================
echo `date` 

TEST_DATA_OUT= /methods/CLAIRE/CLAIRE_code/CLAIRE/Results/claire_test_metanetx.npz
conda activate rxnfp-env
python methods/CLAIRE/CLAIRE_scripts/create_fps.py --query_path "$QUERY_PATH" --fps_path "$FPS_PATH" --test_out "$TEST_DATA_OUT"

echo ===================================
echo ======  Run CLAIRE  ===========
echo ===================================
conda activate claire
python /methods/CLAIRE/CLAIRE_scripts/query_claire.py --test_data_path "$TEST_DATA_OUT" --train_data_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/data/pred_metanetx/lookup_array_ec3.pkl
    / --train_labels_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/data/pred_metanetx/labels_train_ec3.pkl --test_csv_path /data/MetaNetX/test_reactions.csv 
    / --reaction_id_col reaction_id --model_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/results/model/pred_rxn_metanetx/train_final.pth --gmm_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/gmm/gmm_ensumble.pkl

python /methods/CLAIRE/CLAIRE_scripts/get_ec_predictions.py --input_path /methods/CLAIRE/CLAIRE_code/CLAIRE/Results/results_metanetx.csv 
    / --output_path results/Case2/CLAIRE.csv 

################################################### CASE STUDY ###########################################################################

echo ===================================
echo ===  Generate DRFP Fingerprints  ===
echo ===================================
echo `date`

QUERY_PATH=/data/Drugs/reaction_smiles_can.txt
FPS_PATH=/CLAIRE_code/CLAIRE/fps/my_rxn_fps.pkl
conda activate claire
drfp "$QUERY_PATH"  "$FPS_PATH" -d 256
conda deactivate 

echo ================================================
echo ======  Generate fps and concatenate  ===========
echo =================================================
echo `date` 

TEST_DATA_OUT= /methods/CLAIRE/CLAIRE_code/CLAIRE/Results/claire_test_drugs.npz
conda activate rxnfp-env
python methods/CLAIRE/CLAIRE_scripts/create_fps.py --query_path "$QUERY_PATH" --fps_path "$FPS_PATH" --test_out "$TEST_DATA_OUT"

echo ===================================
echo ======  Run CLAIRE  ===========
echo ===================================
conda activate claire
python /methods/CLAIRE/CLAIRE_scripts/query_claire.py --test_data_path "$TEST_DATA_OUT" --train_data_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/data/model_lookup_train.pkl 
    / --train_labels_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/data/pred_rxn_EC123/labels_train_ec3.pkl --test_csv_path /data/Drugs/drug_smiles_updated.csv
    / --reaction_id_col drug --model_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/results/model/pred_rxn_EC123/layer5_node1280_triplet2000_final.pth --gmm_path /methods/CLAIRE/CLAIRE_code/CLAIRE/dev/gmm/gmm_ensumble.pkl

python /methods/CLAIRE/CLAIRE_scripts/get_ec_predictions.py --input_path /methods/CLAIRE/CLAIRE_code/CLAIRE/Results/results_metanetx.csv 
    / --output_path results/CaseStudy/CLAIRE.csv 
