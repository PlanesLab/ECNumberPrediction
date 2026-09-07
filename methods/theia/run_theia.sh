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

module load Miniforge3/24.11.3-2

# NOTE: `conda activate` is unreliable in a non-interactive SLURM batch shell (silently stays on
# the module's base env, e.g. "ModuleNotFoundError: No module named 'pandas'") -- prepend the
# env's bin/ to PATH directly instead, which also makes theia-cli (invoked as a subprocess by
# query_theia.py) resolve correctly.
export PATH="/scratch/jarcagniriv/Envs/theia/bin:$PATH"

python -c 'import sys; print(sys.version_info[:])'
echo "$PATH"

# SEED selects which seed's data to use throughout this script (KEGG-1.8K pool for Case 1,
# MetaNetX split for Case 2); defaults to the original fixed/un-seeded baseline for backward
# compatibility.
SEED="${SEED:-42}"

# --- Case 1: KEGG --- (SKIP_CASE1=1 skips this: Case 1 uses statistical bootstrap of existing
# baseline predictions instead of live reruns; only Case 2 needs a fresh per-seed run)
if [ "${SKIP_CASE1:-0}" != "1" ]; then
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

    python methods/theia/theia_scripts/query_theia.py \
        --query_file "$KEGG_QUERY_FILE" \
        --reaction_ids_file "$KEGG_TEST_CSV" \
        --reaction_id_column "Reaction ID" \
        --output_file "$KEGG_OUT_DIR/theia.csv"
else
    echo "SKIPPING Case 1 (SKIP_CASE1=1)"
fi

# --- Case 2: MetaNetX ---
SKIP_METANETX="${SKIP_METANETX:-0}"
SKIP_CASESTUDY="${SKIP_CASESTUDY:-0}"
if [ "$SKIP_METANETX" != "1" ]; then
if [ -d "data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}" ]; then
    METANETX_TRAIN="data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}/train.tsv"
    METANETX_TEST="data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}/test.tsv"
    METANETX_QUERIES="data/Splits-DBs/MetaNetX/seed_splits/seed${SEED}/queries.txt"
else
    METANETX_TRAIN="data/Splits-DBs/MetaNetX/train.tsv"
    METANETX_TEST="data/Splits-DBs/MetaNetX/test.tsv"
    METANETX_QUERIES="data/Splits-DBs/MetaNetX/queries.txt"
fi
DB_PATH="methods/theia/theia_code/theia/DB_seed${SEED}"
MODEL_PATH="methods/theia/theia_code/theia/models_seed${SEED}"

# Train a custom Theia model on the MetaNetX training split
# encode_split_data.py's main() is invoked via typer.run(main) -- typer only turns parameters
# WITH a default into --flags; path/output_path have no default so they're positional args
# (confirmed from the script's own --help: "Usage: encode_split_data.py [OPTIONS] PATH OUTPUT_PATH").
# --path/--output-path never worked; only --seed (which does have a default) is a real flag.
python methods/theia/theia_code/theia/generateDB/encode_split_data.py \
    "$METANETX_TRAIN" "$DB_PATH" \
    --seed "$SEED"

# encode_split_data.py writes files as a flat "<output-path>-<i>-<ec>-<split>.csv" prefix, not
# into a directory (no os.makedirs / no dir-nesting in its to_csv calls) -- train_all.sh joins
# DATA_PATH and PREFIX with "/", so DATA_PATH must be DB_PATH's *parent* dir, not DB_PATH itself.
bash methods/theia/theia_code/theia/train_all.sh "$(dirname "$DB_PATH")" "$MODEL_PATH" "DB_seed${SEED}-0-ec123" "$SEED"

python methods/theia/theia_scripts/query_theia.py \
    --query_file "$METANETX_QUERIES" \
    --reaction_ids_file "$METANETX_TEST" \
    --reaction_id_column reaction_id \
    --model "$MODEL_PATH/DB_seed${SEED}-0-ec123" \
    --output_file results/Case2/results-DBs/MetaNetX/theia.csv
else
    echo "SKIPPING Case 2 MetaNetX (SKIP_METANETX=1)"
fi

# --- Case 2: Rhea splits (Scaffold/Stratified) ---
# RHEA_SPLIT selects which Rhea split(s) to run (space-separated "Scaffold"/"Stratified");
# defaults to both if unset. Data must already exist under
# data/Splits-Rhea/{Scaffold,Stratified}/seed_splits/seed${SEED}/prepared/{train,test}.tsv,queries.txt
# (written by data/Splits-Rhea/scripts/prepare_seed_split_for_models.py, which renames the raw
# Rhea REACTION_ID/REACTION_SMILES/EC_NUMBER columns to the reaction_id/rxn/ec schema
# encode_split_data.py hardcodes).
RHEA_SPLITS="${RHEA_SPLIT:-Scaffold Stratified}"
for SPLIT in $RHEA_SPLITS; do
    SPLIT_LC=$(echo "$SPLIT" | tr '[:upper:]' '[:lower:]')
    PREPARED_DIR="data/Splits-Rhea/${SPLIT}/seed_splits/seed${SEED}/prepared"
    if [ ! -d "$PREPARED_DIR" ]; then
        echo "SKIPPING Theia Rhea $SPLIT seed $SEED: '$PREPARED_DIR' not found (run prepare_seed_split_for_models.py first)"
        continue
    fi
    RHEA_TRAIN="$PREPARED_DIR/train.tsv"
    RHEA_TEST="$PREPARED_DIR/test.tsv"
    RHEA_QUERIES="$PREPARED_DIR/queries.txt"
    RHEA_DB_PATH="methods/theia/theia_code/theia/DB_rhea_${SPLIT_LC}_seed${SEED}"
    RHEA_MODEL_PATH="methods/theia/theia_code/theia/models_rhea_${SPLIT_LC}_seed${SEED}"
    RHEA_OUT_DIR="results/Case2/results-splits/seed_runs/seed${SEED}/${SPLIT}"
    mkdir -p "$RHEA_OUT_DIR"

    # theia-cli's "source" for this model -- source-split-name = "SOURCE-0-ec123" must equal our
    # trained-model file prefix below (see registration NOTE further down).
    THEIA_SOURCE="DB_rhea_${SPLIT_LC}_seed${SEED}"
    MODEL_PT="$RHEA_MODEL_PATH/${THEIA_SOURCE}-0-ec123.pt"

    # SKIP_THEIA_TRAIN=1 (or the .pt already existing) skips retraining -- training already
    # completed once for this seed/split and is expensive; only the query step needed a fix.
    if [ "${SKIP_THEIA_TRAIN:-0}" != "1" ] && [ ! -f "$MODEL_PT" ]; then
        # see the MetaNetX block's NOTE above -- path/output_path are positional, not --flags
        python methods/theia/theia_code/theia/generateDB/encode_split_data.py \
            "$RHEA_TRAIN" "$RHEA_DB_PATH" \
            --seed "$SEED"

        # see the MetaNetX block's NOTE above -- DATA_PATH must be RHEA_DB_PATH's parent, not itself
        bash methods/theia/theia_code/theia/train_all.sh "$(dirname "$RHEA_DB_PATH")" "$RHEA_MODEL_PATH" "${THEIA_SOURCE}-0-ec123" "$SEED"
    else
        echo "Skipping Theia training for $SPLIT seed $SEED ($MODEL_PT already exists)"
    fi

    # theia-cli's `predict` command loads models from a fixed platformdirs data directory
    # (~/.local/share/theia/) by a "source.name" id, resolved to "<source>-0-<name>.pt" /
    # "-le.pkl" / "-background.pkl" plus a shared "-map.pkl". Symlink the 3 output files
    # into place, plus a stub "-map.pkl" (unused by our code path, but load_models()
    # always unpickles it, so it must exist).
    THEIA_DATA_DIR="$(python3 -c 'import platformdirs; print(platformdirs.user_data_path("theia", "daenuprobst"))')"
    mkdir -p "$THEIA_DATA_DIR"
    for suffix in ".pt" "-le.pkl" "-background.pkl"; do
        ln -sf "$(pwd)/${RHEA_MODEL_PATH}/${THEIA_SOURCE}-0-ec123${suffix}" "$THEIA_DATA_DIR/${THEIA_SOURCE}-0-ec123${suffix}"
    done
    if [ ! -f "$THEIA_DATA_DIR/${THEIA_SOURCE}-map.pkl" ]; then
        python3 -c "import pickle; pickle.dump({}, open('$THEIA_DATA_DIR/${THEIA_SOURCE}-map.pkl', 'wb'))"
    fi

    python methods/theia/theia_scripts/query_theia.py \
        --query_file "$RHEA_QUERIES" \
        --reaction_ids_file "$RHEA_TEST" \
        --reaction_id_column reaction_id \
        --model "${THEIA_SOURCE}.ec123" \
        --output_file "$RHEA_OUT_DIR/theia.csv"
done

# --- Case Study ---
if [ "$SKIP_CASESTUDY" != "1" ]; then
python methods/theia/theia_scripts/query_theia.py \
    --query_file data/Drugs/reaction_smiles_can.txt \
    --reaction_ids_file data/Drugs/drug_smiles_updated.csv \
    --reaction_id_column drug \
    --output_file results/CaseStudy/results/theia.csv
else
    echo "SKIPPING Case Study (SKIP_CASESTUDY=1)"
fi
