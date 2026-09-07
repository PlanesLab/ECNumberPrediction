#!/bin/bash
###################################################################################
### Run BEC-Pred on the Rhea Stratified split (Case 2): baseline (cofactors
### intact) vs. cofactors stripped from BOTH training data and test queries
### (train_nocofactor.tsv/test_nocofactor.tsv, built by
### methods/SIMMER/SIMMER_scripts/strip_cofactors.py). BEC-Pred fine-tunes
### directly on our own train split (unlike SelenzymeRF's external reference
### DB), so this mirrors the SIMMER no-cofactor test almost exactly.
###
### Regenerates baseline fresh rather than reusing results/Case2/results-splits/
### Stratified_becpred.csv -- that file has 1519 rows against a test.tsv that
### only has 757, i.e. it's from a stale/different data vintage and isn't a
### trustworthy comparison point (same kind of mismatch hit earlier this
### session with SIMMER's two evaluation CSVs).
###################################################################################

#SBATCH --job-name=becpred_rhea_nocof
#SBATCH --qos=regular
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --partition=preemption
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jarcagniriv@unav.es
#SBATCH -o results/logs/becpred_rhea_nocof_%j.out

set -euo pipefail

echo "$(date)"
PYTHON="/scratch/jarcagniriv/Envs/becpred_gpu/bin/python3"
$PYTHON -c 'import torch; print(f"GPU: {torch.cuda.is_available()}, count: {torch.cuda.device_count()}")'

CODE_DIR="methods/BEC-Pred/BEC-Pred_code"
SCRIPTS_DIR="methods/BEC-Pred/BEC-Pred_scripts"
RESULTS_DIR="results/Case2/results-splits-nocofactor/BEC-Pred"
mkdir -p "$CODE_DIR/results" results/logs "$RESULTS_DIR"

run_track () {
    local TRACK_NAME=$1     # baseline | nocofactor
    local PREPARED_DIR=$2   # data/Splits-Rhea/Stratified/{baseline,nocofactor}_prepared

    echo "==================================="
    echo "===  BEC-Pred: $TRACK_NAME       ==="
    echo "==================================="
    start_time=$(date +%s)

    local DB_DIR="$CODE_DIR/DB_rhea_stratified_${TRACK_NAME}"
    local MODEL_DIR="$CODE_DIR/model/rhea_stratified_${TRACK_NAME}"

    $PYTHON "$CODE_DIR/generate_becpred_db.py" \
        --train_file "$PREPARED_DIR/train.tsv" \
        --test_file "$PREPARED_DIR/test.tsv" \
        --output_dir "$DB_DIR" \
        --seed 42

    $PYTHON "$CODE_DIR/labels/build_labels_pickle.py" \
        --ec_class_labels_csv "$DB_DIR/ec_class_labels.csv" \
        --output "$DB_DIR/labels_${TRACK_NAME}.pkl"

    $PYTHON "$CODE_DIR/finetune_bec.py" \
        --output_dir "$MODEL_DIR" \
        --train_data "$DB_DIR/train_metanetx.csv" \
        --seed 42

    $PYTHON "$CODE_DIR/eval_model.py" \
        --model_path "$MODEL_DIR" \
        --queries "$PREPARED_DIR/queries.txt" \
        --reaction_ids_file "$PREPARED_DIR/test.tsv" \
        --reaction_id_column reaction_id \
        --output_csv "$CODE_DIR/results/eval_results_rhea_stratified_${TRACK_NAME}.csv"

    $PYTHON "$CODE_DIR/labels/label_assigner.py" \
        --input_csv "$CODE_DIR/results/eval_results_rhea_stratified_${TRACK_NAME}.csv" \
        --labels "$DB_DIR/labels_${TRACK_NAME}.pkl" \
        --output_csv "$RESULTS_DIR/BEC-Pred_${TRACK_NAME}_predictions.csv"

    $PYTHON "$SCRIPTS_DIR/evaluate_rhea_stratified.py" \
        --predictions "$RESULTS_DIR/BEC-Pred_${TRACK_NAME}_predictions.csv" \
        --test-tsv "$PREPARED_DIR/test.tsv" \
        --output "$RESULTS_DIR/BEC-Pred_${TRACK_NAME}_evaluation.csv"

    echo "$TRACK_NAME completed in $(( $(date +%s) - start_time )) seconds"
}

run_track "baseline"   "data/Splits-Rhea/Stratified/baseline_prepared/canonical_becpred"
run_track "nocofactor" "data/Splits-Rhea/Stratified/nocofactor_prepared/canonical_becpred"

echo "==================================="
echo "===  Baseline vs no-cofactor    ==="
echo "==================================="
$PYTHON "$SCRIPTS_DIR/compare_baseline_vs_nocofactor.py" \
    --baseline "$RESULTS_DIR/BEC-Pred_baseline_evaluation.csv" \
    --nocofactor "$RESULTS_DIR/BEC-Pred_nocofactor_evaluation.csv" \
    --output "results/Case2/BECPred_baseline_vs_nocofactor_comparison.csv"

echo "$(date)"
echo "All done."
