#!/bin/bash
###################################################################################
# Author: Josefina Arcagni
# Date: 9/9/2025
# Description: Run BEC-Pred for Case 1, Case 2, and Case Study.
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
#SBATCH -o results/logs/becpred_%j.out

set -euo pipefail

echo "==================================="
echo "===     Load the Packages       ==="
echo "==================================="
echo "$(date)"

module load Anaconda3
conda activate becpred_gpu

python -c 'import sys; print(sys.version_info[:])'
python -c 'import torch; print(f"GPU: {torch.cuda.is_available()}, count: {torch.cuda.device_count()}")'

CODE_DIR="methods/BEC-Pred/BEC-Pred_code"

# --- Case 1: KEGG ---
echo "==================================="
echo "===   Run Case 1 Evaluation     ==="
echo "==================================="
start_time=$(date +%s)

python "$CODE_DIR/eval_model.py" \
    --model_path "$CODE_DIR/model/trained_512" \
    --queries "data/Subsets/KEGG/8:2KEGGTest_canonicalized.txt" \
    --output_csv "$CODE_DIR/results/eval_results_case1.csv"

python "$CODE_DIR/labels/label_assigner.py" \
    --input_csv "$CODE_DIR/results/eval_results_case1.csv" \
    --labels "$CODE_DIR/labels/labels_becpred.pkl" \
    --output_csv "results/Case1/KEGG-1.8K/BEC-Pred.csv"

echo "Case 1 completed in $(( $(date +%s) - start_time )) seconds"

# --- Case 2: MetaNetX ---
echo "==================================="
echo "===   Run Case 2 Evaluation     ==="
echo "==================================="
start_time=$(date +%s)

python "$CODE_DIR/pretrain.py" \
    --output_dir "$CODE_DIR/model/pretrained"

python "$CODE_DIR/finetune_bec.py" \
    --pretrained_model "$CODE_DIR/model/pretrained" \
    --output_dir "$CODE_DIR/model/metanetx" \
    --train_data data/Splits-DBs/MetaNetX/train.tsv

python "$CODE_DIR/eval_model.py" \
    --model_path "$CODE_DIR/model/metanetx" \
    --queries data/Splits-DBs/MetaNetX/queries.txt \
    --output_csv "$CODE_DIR/results/eval_results_case2.csv"

python "$CODE_DIR/labels/label_assigner.py" \
    --input_csv "$CODE_DIR/results/eval_results_case2.csv" \
    --labels "$CODE_DIR/labels/labels_metanetx.pkl" \
    --output_csv "results/Case2/results-DBs/MetaNetX/BEC-Pred.csv"

echo "Case 2 completed in $(( $(date +%s) - start_time )) seconds"

# --- Case Study ---
echo "==================================="
echo "===   Run Case Study Evaluation ==="
echo "==================================="
start_time=$(date +%s)

python "$CODE_DIR/eval_model.py" \
    --model_path "$CODE_DIR/model/trained_512" \
    --queries data/Drugs/reaction_smiles_can.txt \
    --output_csv "$CODE_DIR/results/eval_results_casestudy.csv"

python "$CODE_DIR/labels/label_assigner.py" \
    --input_csv "$CODE_DIR/results/eval_results_casestudy.csv" \
    --labels "$CODE_DIR/labels/labels_becpred.pkl" \
    --output_csv "results/CaseStudy/results/BEC-Pred.csv"

echo "Case Study completed in $(( $(date +%s) - start_time )) seconds"

echo "==================================="
echo "===        Finished Run         ==="
echo "==================================="
echo "$(date)"
