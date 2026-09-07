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

module load Miniforge3/24.11.3-2

# NOTE: `conda activate` is unreliable in a non-interactive SLURM batch shell (silently stays on
# the module's base env) -- invoke the env's python binary directly instead.
PYTHON="/scratch/jarcagniriv/Envs/becpred_gpu/bin/python3"

$PYTHON -c 'import sys; print(sys.version_info[:])'
$PYTHON -c 'import torch; print(f"GPU: {torch.cuda.is_available()}, count: {torch.cuda.device_count()}")'

CODE_DIR="methods/BEC-Pred/BEC-Pred_code"
mkdir -p "$CODE_DIR/results"  # eval_model.py writes here but doesn't create the dir itself

# SEED selects which seed's data to use throughout this script (see run_theia.sh for the same
# convention); defaults to the original fixed/un-seeded baseline for backward compatibility.
SEED="${SEED:-42}"

# SKIP_CASE1/SKIP_METANETX/SKIP_CASESTUDY let a submission scope itself to only the Rhea-splits
# Case 2 block below (e.g. for a Rhea-only bootstrap run) without paying for a redundant KEGG
# eval, MetaNetX retrain, or CaseStudy eval every seed. Same convention as Theia/CLAIRE's
# existing SKIP_CASE1. Default 0 (run everything) for backward compatibility.
SKIP_CASE1="${SKIP_CASE1:-0}"
SKIP_METANETX="${SKIP_METANETX:-0}"
SKIP_CASESTUDY="${SKIP_CASESTUDY:-0}"

# --- Case 1: KEGG ---
if [ "$SKIP_CASE1" != "1" ]; then
echo "==================================="
echo "===   Run Case 1 Evaluation     ==="
echo "==================================="
start_time=$(date +%s)

if [ -d "data/Subsets/KEGG/seed_pools/seed${SEED}" ]; then
    KEGG_QUERY_FILE="data/Subsets/KEGG/seed_pools/seed${SEED}/8:2KEGGTest_canonicalized.txt"
    KEGG_OUT_DIR="results/Case1/seed_runs/seed${SEED}/KEGG-1.8K"
else
    KEGG_QUERY_FILE="data/Subsets/KEGG/8:2KEGGTest_canonicalized.txt"
    KEGG_OUT_DIR="results/Case1/KEGG-1.8K"
fi
mkdir -p "$KEGG_OUT_DIR"

$PYTHON "$CODE_DIR/eval_model.py" \
    --model_path "$CODE_DIR/model/trained_512" \
    --queries "$KEGG_QUERY_FILE" \
    --output_csv "$CODE_DIR/results/eval_results_case1_seed${SEED}.csv"

$PYTHON "$CODE_DIR/labels/label_assigner.py" \
    --input_csv "$CODE_DIR/results/eval_results_case1_seed${SEED}.csv" \
    --labels "$CODE_DIR/labels/labels_becpred.pkl" \
    --output_csv "$KEGG_OUT_DIR/BEC-Pred.csv"

echo "Case 1 completed in $(( $(date +%s) - start_time )) seconds"
else
    echo "SKIPPING Case 1 (SKIP_CASE1=1)"
fi

# --- Case 2: MetaNetX ---
if [ "$SKIP_METANETX" != "1" ]; then
echo "==================================="
echo "===   Run Case 2 Evaluation     ==="
echo "==================================="
start_time=$(date +%s)

if [ -d "data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}" ]; then
    METANETX_TRAIN="data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}/train.tsv"
    METANETX_TEST="data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}/test.tsv"
    METANETX_QUERIES="data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}/queries.txt"
else
    METANETX_TRAIN="data/Splits-DBs/MetaNetX/train.tsv"
    METANETX_TEST="data/Splits-DBs/MetaNetX/test.tsv"
    METANETX_QUERIES="data/Splits-DBs/MetaNetX/queries.txt"
fi
DB_DIR="$CODE_DIR/DB_seed${SEED}"

# NOTE: pretrain.py's MLM pretraining corpus (mlm_train.txt/mlm_test.txt) does not exist
# anywhere in this repo, so pretraining from scratch is not reproducible here. Fine-tune
# directly from rxnfp's shipped bert_pretrained checkpoint instead (finetune_bec.py's
# default --pretrained_model) rather than from model/trained_512 or bert_class_ec_final,
# both of which are already fine-tuned (BertForSequenceClassification) checkpoints, not a
# genuine pretrained-only base model.
$PYTHON "$CODE_DIR/generate_becpred_db.py" \
    --train_file "$METANETX_TRAIN" \
    --test_file "$METANETX_TEST" \
    --output_dir "$DB_DIR" \
    --seed "$SEED"

$PYTHON "$CODE_DIR/labels/build_labels_pickle.py" \
    --ec_class_labels_csv "$DB_DIR/ec_class_labels.csv" \
    --output "$DB_DIR/labels_metanetx_seed${SEED}.pkl"

$PYTHON "$CODE_DIR/finetune_bec.py" \
    --output_dir "$CODE_DIR/model/metanetx_seed${SEED}" \
    --train_data "$DB_DIR/train_metanetx.csv" \
    --seed "$SEED"

$PYTHON "$CODE_DIR/eval_model.py" \
    --model_path "$CODE_DIR/model/metanetx_seed${SEED}" \
    --queries "$METANETX_QUERIES" \
    --output_csv "$CODE_DIR/results/eval_results_case2_seed${SEED}.csv"

$PYTHON "$CODE_DIR/labels/label_assigner.py" \
    --input_csv "$CODE_DIR/results/eval_results_case2_seed${SEED}.csv" \
    --labels "$DB_DIR/labels_metanetx_seed${SEED}.pkl" \
    --output_csv "results/Case2/results-DBs/MetaNetX/BEC-Pred.csv"

echo "Case 2 completed in $(( $(date +%s) - start_time )) seconds"
else
    echo "SKIPPING Case 2 MetaNetX (SKIP_METANETX=1)"
fi

# --- Case 2: Rhea splits (Scaffold/Stratified) ---
# RHEA_SPLIT selects which Rhea split(s) to run (space-separated "Scaffold"/"Stratified");
# defaults to both if unset. Data must already exist under
# data/Splits-Rhea/{Scaffold,Stratified}/seed_splits/seed${SEED}/prepared/{train,test}.tsv,queries.txt
# (written by data/Splits-Rhea/scripts/prepare_seed_split_for_models.py, which renames the raw
# Rhea REACTION_ID/REACTION_SMILES/EC_NUMBER columns to this script's expected reaction_id/
# reaction_smiles/ec schema -- generate_becpred_db.py hardcodes 'reaction_smiles').
echo "==================================="
echo "===   Run Case 2 Evaluation (Rhea splits) ==="
echo "==================================="
start_time=$(date +%s)

RHEA_SPLITS="${RHEA_SPLIT:-Scaffold Stratified}"
for SPLIT in $RHEA_SPLITS; do
    SPLIT_LC=$(echo "$SPLIT" | tr '[:upper:]' '[:lower:]')
    PREPARED_DIR="data/Splits-Rhea/${SPLIT}/seed_splits/seed${SEED}/prepared"
    if [ ! -d "$PREPARED_DIR" ]; then
        echo "SKIPPING Rhea $SPLIT seed $SEED: '$PREPARED_DIR' not found (run prepare_seed_split_for_models.py first)"
        continue
    fi
    # BEC-Pred is a pretrained rxnfp/BERT model -- canonicalize reaction SMILES before
    # training/querying it (same RDKit canonicalization already used for KEGG's
    # *_canonicalized.txt files) so it sees a single consistent SMILES representation
    # instead of whatever arbitrary form Rhea happens to store. canonicalize_for_becpred.py
    # falls back to the original SMILES on any RDKit failure rather than dropping the row,
    # since eval_model.py's queries.txt/test.tsv alignment is purely positional.
    CANON_DIR="$PREPARED_DIR/canonical_becpred"
    # canonicalize_rxn_SMILES.py uses `str | None` type hints (Python 3.10+ syntax) -- $PYTHON
    # here is becpred_gpu's Python 3.6, which can't parse that. Use a working 3.10+ env with
    # rdkit+pandas instead, just for this preprocessing step (unrelated to BEC-Pred's own model).
    /scratch/jarcagniriv/Envs/clean_ec/bin/python3 data/Preprocessing/canonicalize_for_becpred.py \
        --train "$PREPARED_DIR/train.tsv" \
        --test "$PREPARED_DIR/test.tsv" \
        --output_dir "$CANON_DIR" \
        > /dev/null
    RHEA_TRAIN="$CANON_DIR/train.tsv"
    RHEA_TEST="$CANON_DIR/test.tsv"
    RHEA_QUERIES="$CANON_DIR/queries.txt"
    RHEA_DB_DIR="$CODE_DIR/DB_rhea_${SPLIT_LC}_seed${SEED}"
    RHEA_OUT_DIR="results/Case2/results-splits/seed_runs/seed${SEED}/${SPLIT}"
    mkdir -p "$RHEA_OUT_DIR"

    $PYTHON "$CODE_DIR/generate_becpred_db.py" \
        --train_file "$RHEA_TRAIN" \
        --test_file "$RHEA_TEST" \
        --output_dir "$RHEA_DB_DIR" \
        --seed "$SEED"

    $PYTHON "$CODE_DIR/labels/build_labels_pickle.py" \
        --ec_class_labels_csv "$RHEA_DB_DIR/ec_class_labels.csv" \
        --output "$RHEA_DB_DIR/labels_rhea_${SPLIT_LC}_seed${SEED}.pkl"

    # Fine-tune from model/trained_512 (an already EC-classification-adapted BERT checkpoint
    # from an unrelated external project/dataset -- no leakage risk against our Rhea test sets)
    # instead of rxnfp's generic reaction-pretrained-only checkpoint, as a warm start for
    # transfer learning. trained_512's own label count (308) differs from our num_labels=353,
    # so its classification head gets reinitialized fresh on load while the BERT encoder
    # weights (already adapted toward EC-relevant reaction representations) carry over.
    $PYTHON "$CODE_DIR/finetune_bec.py" \
        --output_dir "$CODE_DIR/model/rhea_${SPLIT_LC}_seed${SEED}" \
        --train_data "$RHEA_DB_DIR/train_metanetx.csv" \
        --pretrained_model "$CODE_DIR/model/trained_512" \
        --seed "$SEED"

    # --reaction_ids_file/--reaction_id_column attach a 'reaction_id' column to the output
    # (eval_model.py otherwise writes only a bare Prediction column with no ID at all, which
    # join_results.py cannot merge against ground truth) -- RHEA_TEST is row-order-aligned with
    # RHEA_QUERIES per prepare_seed_split_for_models.py's contract.
    $PYTHON "$CODE_DIR/eval_model.py" \
        --model_path "$CODE_DIR/model/rhea_${SPLIT_LC}_seed${SEED}" \
        --queries "$RHEA_QUERIES" \
        --reaction_ids_file "$RHEA_TEST" \
        --reaction_id_column reaction_id \
        --output_csv "$CODE_DIR/results/eval_results_case2_rhea_${SPLIT_LC}_seed${SEED}.csv"

    $PYTHON "$CODE_DIR/labels/label_assigner.py" \
        --input_csv "$CODE_DIR/results/eval_results_case2_rhea_${SPLIT_LC}_seed${SEED}.csv" \
        --labels "$RHEA_DB_DIR/labels_rhea_${SPLIT_LC}_seed${SEED}.pkl" \
        --output_csv "$RHEA_OUT_DIR/BEC-Pred.csv"
done

echo "Case 2 (Rhea splits) completed in $(( $(date +%s) - start_time )) seconds"

# --- Case Study ---
if [ "$SKIP_CASESTUDY" != "1" ]; then
echo "==================================="
echo "===   Run Case Study Evaluation ==="
echo "==================================="
start_time=$(date +%s)

$PYTHON "$CODE_DIR/eval_model.py" \
    --model_path "$CODE_DIR/model/trained_512" \
    --queries data/Drugs/reaction_smiles_can.txt \
    --output_csv "$CODE_DIR/results/eval_results_casestudy.csv"

$PYTHON "$CODE_DIR/labels/label_assigner.py" \
    --input_csv "$CODE_DIR/results/eval_results_casestudy.csv" \
    --labels "$CODE_DIR/labels/labels_becpred.pkl" \
    --output_csv "results/CaseStudy/results/BEC-Pred.csv"

echo "Case Study completed in $(( $(date +%s) - start_time )) seconds"
else
    echo "SKIPPING Case Study (SKIP_CASESTUDY=1)"
fi

echo "==================================="
echo "===        Finished Run         ==="
echo "==================================="
echo "$(date)"
