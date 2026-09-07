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

module load Miniforge3/24.11.3-2
python -c 'import sys; print(sys.version_info[:])'
echo "$PATH"

CLAIRE_CODE="methods/CLAIRE/CLAIRE_code/CLAIRE"
CLAIRE_SCRIPTS="methods/CLAIRE/CLAIRE_scripts"

# drfp writes its output pickle without creating the parent dir first (silently no-ops if
# fps/ is missing, which only failed to surface until a seed other than the original 42 was
# used, since fps/ happened to already exist from that earlier run).
mkdir -p "$CLAIRE_CODE/fps"

# --- Case 1: KEGG --- (SKIP_CASE1=1 skips this: Case 1 uses statistical bootstrap of existing
# baseline predictions instead of live reruns; only Case 2 needs a fresh per-seed run)
# NOTE: requires dev/data/model_lookup_train.pkl + dev/data/pred_rxn_EC123/labels_train_ec3.pkl,
# which upstream CLAIRE distributes separately via Zenodo (https://zenodo.org/records/14635841),
# not through the git clone -- download and place under $CLAIRE_CODE/dev/ before running this block.
SEED="${SEED:-42}"
if [ "${SKIP_CASE1:-0}" != "1" ]; then
    echo "==================================="
    echo "===  Generate DRFP Fingerprints  ==="
    echo "==================================="
    echo "$(date)"

    if [ -d "data/Subsets/KEGG/seed_pools/seed${SEED}" ]; then
        KEGG_QUERY_FILE="data/Subsets/KEGG/seed_pools/seed${SEED}/8:2KEGGTest_canonicalized.txt"
        KEGG_TEST_CSV="data/Subsets/KEGG/seed_pools/seed${SEED}/kegg_reactions_current_test.csv"
        KEGG_OUT_DIR="results/Case1/seed_runs/seed${SEED}/KEGG-1.8K"
    else
        KEGG_QUERY_FILE="data/Subsets/KEGG/8:2KEGGTest_canonicalized.txt"
        KEGG_TEST_CSV="data/Subsets/KEGG/kegg_reactions_current_test.csv"
        KEGG_OUT_DIR="results/Case1/KEGG-1.8K"
    fi
    mkdir -p "$KEGG_OUT_DIR"

    QUERY_PATH="$KEGG_QUERY_FILE"
    FPS_PATH="$CLAIRE_CODE/fps/my_rxn_fps_kegg_seed${SEED}.pkl"
    TEST_DATA_OUT="$CLAIRE_CODE/Results/claire_test_kegg_seed${SEED}.npz"

    export PATH="/scratch/jarcagniriv/Envs/claire/bin:$PATH"
    /scratch/jarcagniriv/Envs/theia/bin/drfp "$QUERY_PATH" "$FPS_PATH" -d 256

    export PATH="/scratch/jarcagniriv/Envs/rxnfp-env/bin:$PATH"
    export LD_LIBRARY_PATH="/scratch/jarcagniriv/Envs/rxnfp-env/lib:${LD_LIBRARY_PATH:-}"
    /scratch/jarcagniriv/Envs/rxnfp-env/bin/python "$CLAIRE_SCRIPTS/create_fps.py" \
        "$QUERY_PATH" "$FPS_PATH" "$TEST_DATA_OUT"

    export PATH="/scratch/jarcagniriv/Envs/claire/bin:$PATH"
    /scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_SCRIPTS/query_claire.py" \
        --test_data_path "$TEST_DATA_OUT" \
        --train_data_path "$CLAIRE_CODE/dev/data/model_lookup_train.pkl" \
        --train_labels_path "$CLAIRE_CODE/dev/data/pred_rxn_EC123/labels_train_ec3.pkl" \
        --test_csv_path "$KEGG_TEST_CSV" \
        --reaction_id_col "Reaction ID" \
        --model_path "$CLAIRE_CODE/dev/results/model/pred_rxn_EC123/layer5_node1280_triplet2000_final.pth" \
        --gmm_path "$CLAIRE_CODE/dev/gmm/gmm_ensumble.pkl" \
        --out_filename "$CLAIRE_CODE/Results/results_kegg_seed${SEED}"

    /scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_SCRIPTS/get_ec_predictions.py" \
        --input "$CLAIRE_CODE/Results/results_kegg_seed${SEED}_prediction.csv" \
        --output "$KEGG_OUT_DIR/CLAIRE.csv"
else
    echo "SKIPPING Case 1 (SKIP_CASE1=1)"
fi

# --- Case 2: MetaNetX ---
# SEED selects which seed's train/test split to use (see run_theia.sh/run_becpred.sh for the same
# convention); defaults to the original un-seeded baseline split for backward compatibility.
# Case 2's embedding DB (esm_emb_dict_ec3.pkl/lookup_array_ec3.pkl/labels_train_ec3.pkl) is built
# entirely from MetaNetX data below -- no Zenodo dependency, unlike Case 1.
SKIP_METANETX="${SKIP_METANETX:-0}"
SKIP_CASESTUDY="${SKIP_CASESTUDY:-0}"
SEED="${SEED:-42}"
# EPOCHS is set here (not inside the SKIP_METANETX guard below) because the Rhea-splits block
# further down also needs it, and would hit an unbound-variable error under set -u if
# SKIP_METANETX=1 skipped this assignment.
EPOCHS="${CLAIRE_EPOCHS:-2000}"  # upstream's own default. The 200-epoch cut previously used
# here was found to cause severe seed-to-seed instability (one seed/split combo scored MCC=0.003,
# essentially random, while others on the same split scored 0.24-0.27) -- root cause: CLAIRE's
# adaptive_rate=200 hard-negative re-mining reset (upstream default, see
# dev/training/train-pred_rxn_EC.py) never got to complete even a single full cycle at only 200
# epochs. Verified: rerunning the MCC=0.003 seed/split at 2000 epochs raised it to 0.298, in line
# with its peers. Each epoch is only ~0.1-0.15s, so the extra cost is a few minutes, not hours.
if [ "$SKIP_METANETX" != "1" ]; then
echo "==================================="
echo "===  Train CLAIRE model  ==="
echo "==================================="
echo "$(date)"

if [ -d "data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}" ]; then
    METANETX_TRAIN="data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}/train.tsv"
    METANETX_TEST="data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}/test.tsv"
    METANETX_QUERIES="data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}/queries.txt"
else
    METANETX_TRAIN="data/Splits-DBs/MetaNetX/train.tsv"
    METANETX_TEST="data/Splits-DBs/MetaNetX/test.tsv"
    METANETX_QUERIES="data/Splits-DBs/MetaNetX/queries.txt"
fi
DB_DIR="$CLAIRE_CODE/dev/data/pred_metanetx_seed${SEED}"
MODEL_DIR="$CLAIRE_CODE/dev/results/model/pred_rxn_metanetx_seed${SEED}"

export PATH="/scratch/jarcagniriv/Envs/rxnfp-env/bin:$PATH"
export LD_LIBRARY_PATH="/scratch/jarcagniriv/Envs/rxnfp-env/lib:${LD_LIBRARY_PATH:-}"
/scratch/jarcagniriv/Envs/rxnfp-env/bin/python "$CLAIRE_CODE/rxnfp_create.py" \
    --train_file "$METANETX_TRAIN" \
    --output_dir "$DB_DIR"

export PATH="/scratch/jarcagniriv/Envs/claire/bin:$PATH"
/scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_CODE/create_claire_db.py" \
    --db_dir "$DB_DIR"

/scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_CODE/dev/training/train-pred_rxn_EC.py" \
    --db_dir "$DB_DIR" \
    --output_dir "$MODEL_DIR" \
    --epoch "$EPOCHS" \
    --seed "$SEED"

QUERY_PATH="$METANETX_QUERIES"
FPS_PATH="$CLAIRE_CODE/fps/my_rxn_fps_seed${SEED}.pkl"
TEST_DATA_OUT="$CLAIRE_CODE/Results/claire_test_metanetx_seed${SEED}.npz"

/scratch/jarcagniriv/Envs/theia/bin/drfp "$QUERY_PATH" "$FPS_PATH" -d 256

export PATH="/scratch/jarcagniriv/Envs/rxnfp-env/bin:$PATH"
export LD_LIBRARY_PATH="/scratch/jarcagniriv/Envs/rxnfp-env/lib:${LD_LIBRARY_PATH:-}"
/scratch/jarcagniriv/Envs/rxnfp-env/bin/python "$CLAIRE_SCRIPTS/create_fps.py" \
    "$QUERY_PATH" "$FPS_PATH" "$TEST_DATA_OUT"

export PATH="/scratch/jarcagniriv/Envs/claire/bin:$PATH"
/scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_SCRIPTS/query_claire.py" \
    --test_data_path "$TEST_DATA_OUT" \
    --train_data_path "$DB_DIR/lookup_array_ec3.pkl" \
    --train_labels_path "$DB_DIR/labels_train_ec3.pkl" \
    --test_csv_path "$METANETX_TEST" \
    --reaction_id_col reaction_id \
    --model_path "$MODEL_DIR/train_final.pth" \
    --gmm_path "$CLAIRE_CODE/dev/gmm/gmm_ensumble.pkl" \
    --out_filename "$CLAIRE_CODE/Results/results_metanetx_seed${SEED}"

/scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_SCRIPTS/get_ec_predictions.py" \
    --input "$CLAIRE_CODE/Results/results_metanetx_seed${SEED}_prediction.csv" \
    --output results/Case2/results-DBs/MetaNetX/CLAIRE.csv
else
    echo "SKIPPING Case 2 MetaNetX (SKIP_METANETX=1)"
fi

# --- Case 2: Rhea splits (Scaffold/Stratified) ---
# RHEA_SPLIT selects which Rhea split(s) to run (space-separated "Scaffold"/"Stratified");
# defaults to both if unset. Data must already exist under
# data/Splits-Rhea/{Scaffold,Stratified}/seed_splits/seed${SEED}/prepared/{train,test}.tsv,queries.txt
# (written by data/Splits-Rhea/scripts/prepare_seed_split_for_models.py, which renames the raw
# Rhea REACTION_ID/REACTION_SMILES/EC_NUMBER columns to reaction_id/reaction_smiles/ec -- the
# defaults rxnfp_create.py/create_claire_db.py already expect).
echo "==================================="
echo "===  Train CLAIRE model (Rhea splits)  ==="
echo "==================================="
echo "$(date)"

RHEA_SPLITS="${RHEA_SPLIT:-Scaffold Stratified}"
for SPLIT in $RHEA_SPLITS; do
    SPLIT_LC=$(echo "$SPLIT" | tr '[:upper:]' '[:lower:]')
    PREPARED_DIR="data/Splits-Rhea/${SPLIT}/seed_splits/seed${SEED}/prepared"
    if [ ! -d "$PREPARED_DIR" ]; then
        echo "SKIPPING CLAIRE Rhea $SPLIT seed $SEED: '$PREPARED_DIR' not found (run prepare_seed_split_for_models.py first)"
        continue
    fi
    RHEA_TRAIN="$PREPARED_DIR/train.tsv"
    RHEA_TEST="$PREPARED_DIR/test.tsv"
    RHEA_QUERIES="$PREPARED_DIR/queries.txt"
    RHEA_DB_DIR="$CLAIRE_CODE/dev/data/pred_rhea_${SPLIT_LC}_seed${SEED}"
    RHEA_MODEL_DIR="$CLAIRE_CODE/dev/results/model/pred_rxn_rhea_${SPLIT_LC}_seed${SEED}"
    RHEA_OUT_DIR="results/Case2/results-splits/seed_runs/seed${SEED}/${SPLIT}"
    mkdir -p "$RHEA_OUT_DIR"

    export PATH="/scratch/jarcagniriv/Envs/rxnfp-env/bin:$PATH"
    export LD_LIBRARY_PATH="/scratch/jarcagniriv/Envs/rxnfp-env/lib:${LD_LIBRARY_PATH:-}"
    /scratch/jarcagniriv/Envs/rxnfp-env/bin/python "$CLAIRE_CODE/rxnfp_create.py" \
        --train_file "$RHEA_TRAIN" \
        --output_dir "$RHEA_DB_DIR"

    export PATH="/scratch/jarcagniriv/Envs/claire/bin:$PATH"
    /scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_CODE/create_claire_db.py" \
        --db_dir "$RHEA_DB_DIR"

    /scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_CODE/dev/training/train-pred_rxn_EC.py" \
        --db_dir "$RHEA_DB_DIR" \
        --output_dir "$RHEA_MODEL_DIR" \
        --epoch "$EPOCHS" \
        --seed "$SEED"

    # train-pred_rxn_EC.py and inference_EC.py define separate LayerNormNet classes with
    # different state_dict key naming, so query_claire.py fails to load a raw training
    # checkpoint. Remap before use.
    /scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_SCRIPTS/remap_checkpoint.py" \
        --input "$RHEA_MODEL_DIR/train_final.pth" \
        --output "$RHEA_MODEL_DIR/train_final_remapped.pth"

    QUERY_PATH="$RHEA_QUERIES"
    FPS_PATH="$CLAIRE_CODE/fps/my_rxn_fps_rhea_${SPLIT_LC}_seed${SEED}.pkl"
    TEST_DATA_OUT="$CLAIRE_CODE/Results/claire_test_rhea_${SPLIT_LC}_seed${SEED}.npz"

    /scratch/jarcagniriv/Envs/theia/bin/drfp "$QUERY_PATH" "$FPS_PATH" -d 256

    export PATH="/scratch/jarcagniriv/Envs/rxnfp-env/bin:$PATH"
    export LD_LIBRARY_PATH="/scratch/jarcagniriv/Envs/rxnfp-env/lib:${LD_LIBRARY_PATH:-}"
    /scratch/jarcagniriv/Envs/rxnfp-env/bin/python "$CLAIRE_SCRIPTS/create_fps.py" \
        "$QUERY_PATH" "$FPS_PATH" "$TEST_DATA_OUT"

    export PATH="/scratch/jarcagniriv/Envs/claire/bin:$PATH"
    /scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_SCRIPTS/query_claire.py" \
        --test_data_path "$TEST_DATA_OUT" \
        --train_data_path "$RHEA_DB_DIR/lookup_array_ec3.pkl" \
        --train_labels_path "$RHEA_DB_DIR/labels_train_ec3.pkl" \
        --test_csv_path "$RHEA_TEST" \
        --reaction_id_col reaction_id \
        --model_path "$RHEA_MODEL_DIR/train_final_remapped.pth" \
        --gmm_path "$CLAIRE_CODE/dev/gmm/gmm_ensumble.pkl" \
        --out_filename "$CLAIRE_CODE/Results/results_rhea_${SPLIT_LC}_seed${SEED}"

    /scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_SCRIPTS/get_ec_predictions.py" \
        --input "$CLAIRE_CODE/Results/results_rhea_${SPLIT_LC}_seed${SEED}_prediction.csv" \
        --output "$RHEA_OUT_DIR/CLAIRE.csv"
done

# --- Case Study ---
if [ "$SKIP_CASESTUDY" != "1" ]; then
QUERY_PATH="data/Drugs/reaction_smiles_can.txt"
FPS_PATH="$CLAIRE_CODE/fps/my_rxn_fps.pkl"
TEST_DATA_OUT="$CLAIRE_CODE/Results/claire_test_drugs.npz"

/scratch/jarcagniriv/Envs/theia/bin/drfp "$QUERY_PATH" "$FPS_PATH" -d 256

export PATH="/scratch/jarcagniriv/Envs/rxnfp-env/bin:$PATH"
export LD_LIBRARY_PATH="/scratch/jarcagniriv/Envs/rxnfp-env/lib:${LD_LIBRARY_PATH:-}"
/scratch/jarcagniriv/Envs/rxnfp-env/bin/python "$CLAIRE_SCRIPTS/create_fps.py" \
    "$QUERY_PATH" "$FPS_PATH" "$TEST_DATA_OUT"

export PATH="/scratch/jarcagniriv/Envs/claire/bin:$PATH"
/scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_SCRIPTS/query_claire.py" \
    --test_data_path "$TEST_DATA_OUT" \
    --train_data_path "$CLAIRE_CODE/dev/data/model_lookup_train.pkl" \
    --train_labels_path "$CLAIRE_CODE/dev/data/pred_rxn_EC123/labels_train_ec3.pkl" \
    --test_csv_path data/Drugs/drug_smiles_updated.csv \
    --reaction_id_col drug \
    --model_path "$CLAIRE_CODE/dev/results/model/pred_rxn_EC123/layer5_node1280_triplet2000_final.pth" \
    --gmm_path "$CLAIRE_CODE/dev/gmm/gmm_ensumble.pkl" \
    --out_filename "$CLAIRE_CODE/Results/results_drugs"

/scratch/jarcagniriv/Envs/claire/bin/python "$CLAIRE_SCRIPTS/get_ec_predictions.py" \
    --input "$CLAIRE_CODE/Results/results_drugs_prediction.csv" \
    --output results/CaseStudy/results/CLAIRE.csv
else
    echo "SKIPPING Case Study (SKIP_CASESTUDY=1)"
fi
