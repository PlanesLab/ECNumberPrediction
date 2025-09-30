#!/bin/bash

###################################################################################
# Author: Josefina Arcagni
# Date: 9/9/2025
# Description: Run BEC-Pred for various cases.
###################################################################################

###################################################################################
### Run BEC-Pred 
###################################################################################

#SBATCH --job-name=Becpred
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jarcagniriv@unav.es
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --partition=preemption
#SBATCH -o /scratch/jarcagniriv/ECNumberPrediction/results/logs/bench_%j.out

echo ===================================
echo ===     Load the Packages       ===
echo ===================================
echo `date`

module load Anaconda3  
source activate /scratch/jarcagniriv/ECNumberPredictionReview/Envs/becpred_gpu

python -c 'import sys; print(sys.version_info[:])'

echo $PATH
export PATH="ECNumberPredictionReview/Envs/becpred_gpu/bin:$PATH"
export LD_LIBRARY_PATH="/ECNumberPredictionReview/Envs/becpred_gpu/lib:$LD_LIBRARY_PATH"

# Check GPU availability
python -c 'import torch; print(f"GPU Available: {torch.cuda.is_available()}"); print(f"GPU Count: {torch.cuda.device_count()}"); print(f"GPU Name: {torch.cuda.get_device_name(0)}")'

###################################################################################
### CASE 1
###################################################################################
echo ===================================
echo ===   Run Case1 Evaluation     ===
echo ===================================
start_time=$(date +%s)

CODE_DIR="ECNumberPrediction/methods/BEC-Pred/BEC-Pred_code"
RESULTS_DIR="ECNumberPrediction/results/Case1"
DATA_DIR="ECNumberPrediction/data/KEGG"

EVAL_OUTPUT="$CODE_DIR/results/eval_results.csv"
FINETUNE_OUTPUT="$CODE_DIR/model/trained_512"
QUERIES="$DATA_DIR/8:2KEGGTest_canonicalized.txt"

python $CODE_DIR/eval_model.py --model_path "$FINETUNE_OUTPUT" --queries "$QUERIES" --output_csv "$EVAL_OUTPUT"

end_time=$(date +%s)
eval_duration=$((end_time - start_time))
echo "Case 1 Eval completed in $eval_duration seconds"

echo "Assigning labels..."
FINAL_RESULTS="$RESULTS_DIR/BEC-Pred.csv"
python $CODE_DIR/labels/label_assigner.py --input_csv "$EVAL_OUTPUT" --labels CODE_DIR/labels/labels_becpred.pkl --output_csv "$FINAL_RESULTS"

###################################################################################
### CASE 2
###################################################################################
echo ===================================
echo ===   Run Case2 Evaluation     ===
echo ===================================
start_time=$(date +%s)

CODE_DIR="ECNumberPrediction/methods/BEC-Pred/BEC-Pred_code"
RESULTS_DIR="ECNumberPrediction/results/Case2"
DATA_DIR="ECNumberPrediction/data/MetaNetX"

start_time=$(date +%s)  
python $CODE_DIR/pretrain.py --output_dir $CODE_DIR/model/pretrained
end_time=$(date +%s)  
pretrain_duration=$((end_time - start_time)) 
echo "Pretrain completed in $pretrain_duration seconds"

start_time=$(date +%s)  
python $CODE_DIR/finetune_bec.py --pretrained_model $CODE_DIR/model/pretrained --output_dir $CODE_DIR/model/metanetx --train_data $DATA_DIR/train_metanetx.csv
end_time=$(date +%s)  
finetune_duration=$((end_time - start_time))  
echo "Finetune completed in $finetune_duration seconds"

EVAL_OUTPUT="$CODE_DIR/results/eval_results.csv"
FINETUNE_OUTPUT="$CODE_DIR/model/metanetx"
QUERIES="$DATA_DIR/queries.txt"

python $CODE_DIR/eval_model.py --model_path "$FINETUNE_OUTPUT" --queries "$QUERIES" --output_csv "$EVAL_OUTPUT"

end_time=$(date +%s)
eval_duration=$((end_time - start_time))
echo "Case 2 Eval completed in $eval_duration seconds"

echo "Assigning labels..."
FINAL_RESULTS="$RESULTS_DIR/BEC-Pred.csv"
python $CODE_DIR/labels/label_assigner.py --input_csv "$EVAL_OUTPUT" --labels CODE_DIR/labels/labels_metanetx.pkl --output_csv "$FINAL_RESULTS"

###################################################################################
### CASE STUDY
###################################################################################
echo ===================================
echo ===   Run Case Study Evaluation ===
echo ===================================
start_time=$(date +%s)

CODE_DIR="ECNumberPrediction/methods/BEC-Pred/BEC-Pred_code"
RESULTS_DIR="ECNumberPrediction/results/CaseStudy"
DATA_DIR="ECNumberPrediction/data/Drugs"

EVAL_OUTPUT="$CODE_DIR/results/eval_results.csv"
FINETUNE_OUTPUT="$CODE_DIR/model/trained_512"
QUERIES="$DATA_DIR/reaction_smiles_can.txt"

python $CODE_DIR/eval_model.py --model_path "$FINETUNE_OUTPUT" --queries "$QUERIES" --output_csv "$EVAL_OUTPUT"

end_time=$(date +%s)
eval_duration=$((end_time - start_time))
echo "Case Study Eval completed in $eval_duration seconds"

echo "Assigning labels..."
FINAL_RESULTS="$RESULTS_DIR/BEC-Pred.csv"
python $CODE_DIR/labels/label_assigner.py --input_csv "$EVAL_OUTPUT" --labels CODE_DIR/labels/labels_becpred.pkl --output_csv "$FINAL_RESULTS"

echo ===================================
echo ===        Finished Run         ===
echo ===================================
echo `date`
