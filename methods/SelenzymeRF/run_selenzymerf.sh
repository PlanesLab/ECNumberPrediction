#!/bin/bash
###################################################################################
# Author: Josefina Arcagni
# Date: 9/9/2025
# Description: Run SelenzymeRF for Case 1, Case 2, and Case Study.
###################################################################################
set -euo pipefail

# SKIP_CASE1/SKIP_CASESTUDY let a submission scope itself to only Case 2 (each of Case 1/Case
# Study costs a full ~5-minute Apptainer server start on top of query time, so skipping them
# matters more here than for the other methods). Default 0 (run everything) for compatibility.
SEED="${SEED:-42}"
SKIP_CASE1="${SKIP_CASE1:-0}"
SKIP_CASESTUDY="${SKIP_CASESTUDY:-0}"

# --- Case 1: KEGG ---
# SEED selects which seed's resampled KEGG-1.8K test set to use (see run_theia.sh for the same
# convention); defaults to the original fixed test set for backward compatibility.
if [ "$SKIP_CASE1" != "1" ]; then
echo "Running SelenzymeRF for Case 1"
if [ -d "data/Subsets/KEGG/seed_pools/seed${SEED}" ]; then
    KEGG_TEST_CSV="data/Subsets/KEGG/seed_pools/seed${SEED}/kegg_reactions_current_test.csv"
    KEGG_OUT_DIR="results/Case1/seed_runs/seed${SEED}/KEGG-1.8K"
else
    KEGG_TEST_CSV="data/Subsets/KEGG/kegg_reactions_current_test.csv"
    KEGG_OUT_DIR="results/Case1/KEGG-1.8K"
fi
mkdir -p "$KEGG_OUT_DIR"

# start_server.sh is interactive and would hang in a SLURM batch job; start_server_noninteractive.sh
# runs the same workflow via env vars. SELENZYME2_DIR is seed-specific so parallel bootstrap
# seeds don't collide on the same mount dir (see the Case 2 block below).
SELENZYME2_DIR="selenzyme2_case1_seed${SEED}" \
    bash methods/SelenzymeRF/SelenzymeRF_code/start_server_noninteractive.sh &
SERVER_PID=$!

# Wait for server to start
sleep 300

/scratch/jarcagniriv/Envs/clean_ec/bin/python3 methods/SelenzymeRF/SelenzymeRF_scripts/query_selenzyme.py \
    --csv_file "$KEGG_TEST_CSV" \
    --results_folder "results/SelenzymeRF/case1_seed${SEED}" \
    --server_url http://localhost:5001 \
    --reaction_name_column "Reaction ID" \
    --reaction_smiles_column Isomeric_SMILES

/scratch/jarcagniriv/Envs/clean_ec/bin/python3 methods/SelenzymeRF/SelenzymeRF_scripts/get_results.py \
    --folder "results/SelenzymeRF/case1_seed${SEED}" \
    --ec_column "EC Number" \
    --output "$KEGG_OUT_DIR/SelenzymeRF.csv"

kill "$SERVER_PID" 2>/dev/null || true
else
    echo "SKIPPING Case 1 (SKIP_CASE1=1)"
fi

# --- Case Study ---
if [ "$SKIP_CASESTUDY" != "1" ]; then
echo "Running SelenzymeRF for Case Study"
SELENZYME2_DIR="selenzyme2_casestudy" \
    bash methods/SelenzymeRF/SelenzymeRF_code/start_server_noninteractive.sh &
SERVER_PID=$!

sleep 300

/scratch/jarcagniriv/Envs/clean_ec/bin/python3 methods/SelenzymeRF/SelenzymeRF_scripts/query_selenzyme.py \
    --csv_file data/Drugs/drug_smiles_updated.csv \
    --results_folder results/SelenzymeRF/casestudy \
    --server_url http://localhost:5001 \
    --reaction_name_column drug \
    --reaction_smiles_column reaction_smiles

/scratch/jarcagniriv/Envs/clean_ec/bin/python3 methods/SelenzymeRF/SelenzymeRF_scripts/get_results.py \
    --folder results/SelenzymeRF/casestudy \
    --ec_column "EC Number" \
    --output results/CaseStudy/results/SelenzymeRF.csv

kill "$SERVER_PID" 2>/dev/null || true
else
    echo "SKIPPING Case Study (SKIP_CASESTUDY=1)"
fi

# --- Case 2: Rhea-Scaffold / Rhea-Stratified splits, with a per-seed filtered reference DB ---
# (Rhea-Time is intentionally excluded: it is not reseedable, see data/Splits-Rhea/Time/.)
# (MetaNetX Case 2 is intentionally out of scope here -- generate_selenzyme_db.py supports
#  --id_scheme metanetx too, but building/running that track is left for a later pass.)
echo "Running SelenzymeRF for Case 2 (Rhea Scaffold/Stratified, seed${SEED})"

SELENZYME_CODE_DIR="methods/SelenzymeRF/SelenzymeRF_code"

# Extract the base (unfiltered) MetaNetX reference DB once and reuse it for every
# seed/track -- generate_selenzyme_db.py filters a copy of it per seed, it never modifies
# this directory. Reuses the exact directory name (data_2023/) already listed in
# SelenzymeRF_code/.gitignore.
BASE_DB_DIR="$SELENZYME_CODE_DIR/data_2023"
if [ ! -f "$BASE_DB_DIR/reac_prop.tsv" ]; then
    echo "Extracting base reference DB to $BASE_DB_DIR"
    mkdir -p "$BASE_DB_DIR"
    unzip -q "$SELENZYME_CODE_DIR/compressed_data/data_2023.zip" -d "$BASE_DB_DIR"
    mv "$BASE_DB_DIR"/data_2023/* "$BASE_DB_DIR"/
    rmdir "$BASE_DB_DIR"/data_2023
fi
if [ ! -f "$BASE_DB_DIR/seqs.fasta" ]; then
    echo "Extracting seqs.fasta to $BASE_DB_DIR"
    unzip -q "$SELENZYME_CODE_DIR/compressed_data/seqs.zip" -d "$BASE_DB_DIR"
fi
SEQS_FASTA="$BASE_DB_DIR/seqs.fasta"

RHEA_TRACKS_TO_RUN="${RHEA_SPLIT:-Scaffold Stratified}"  # override to e.g. "Stratified" to rerun just one track
for RHEA_TRACK in $RHEA_TRACKS_TO_RUN; do
    RHEA_TEST_TSV="data/Splits-Rhea/${RHEA_TRACK}/seed_splits/seed${SEED}/test.tsv"
    RHEA_TRAIN_TSV="data/Splits-Rhea/${RHEA_TRACK}/seed_splits/seed${SEED}/train.tsv"
    if [ ! -f "$RHEA_TEST_TSV" ]; then
        echo "SKIPPING SelenzymeRF Case 2 (${RHEA_TRACK}, seed${SEED}): $RHEA_TEST_TSV not found"
        continue
    fi

    # DB_BUILD_MODE=train_only (new default) builds the reference DB from ONLY the seed's train
    # split (--filter_mode include_only against train.tsv), methodologically consistent with how
    # SIMMER's reference DB is built purely from train -- rather than the original approach of
    # excluding just the test reactions from the full external MetaNetX reference (still leakage-
    # safe, but a much broader, not train-scoped, reference pool). Set DB_BUILD_MODE=exclude_test
    # to fall back to the original broad-pool behavior.
    DB_BUILD_MODE="${DB_BUILD_MODE:-train_only}"
    if [ "$DB_BUILD_MODE" == "train_only" ]; then
        FILTER_MODE="include_only"
        FILTER_SOURCE_TSV="$RHEA_TRAIN_TSV"
        FILTERED_DB_DIR="$SELENZYME_CODE_DIR/filtered_dbs/Rhea_${RHEA_TRACK}_seed${SEED}_trainonly"
    else
        FILTER_MODE="exclude"
        FILTER_SOURCE_TSV="$RHEA_TEST_TSV"
        FILTERED_DB_DIR="$SELENZYME_CODE_DIR/filtered_dbs/Rhea_${RHEA_TRACK}_seed${SEED}"
    fi
    echo "Building ${DB_BUILD_MODE} reference DB for Rhea ${RHEA_TRACK} seed${SEED} -> $FILTERED_DB_DIR"
    /scratch/jarcagniriv/Envs/clean_ec/bin/python3 methods/SelenzymeRF/SelenzymeRF_scripts/generate_selenzyme_db.py \
        --base_db_dir "$BASE_DB_DIR" \
        --test_tsv "$FILTER_SOURCE_TSV" \
        --filter_mode "$FILTER_MODE" \
        --id_scheme rhea \
        --output_dir "$FILTERED_DB_DIR" \
        --seqs_fasta "$SEQS_FASTA"

    # start_server_noninteractive.sh cd's into its own script directory before consulting
    # DATA_SOURCE_DIR, so a repo-root-relative path like $FILTERED_DB_DIR silently resolves to
    # a doubled, nonexistent path there ("does not exist" despite generate_selenzyme_db.py having
    # just written it) -- pass an absolute path instead.
    # SELENZYME2_DIR is seed-specific: the default "selenzyme2" mount dir lives under the shared
    # SelenzymeRF_code/ filesystem path, so concurrent seed jobs (this script is submitted 3x in
    # parallel, once per seed) collide there (confirmed: "mkdir: cannot create directory
    # 'selenzyme2': File exists" from a sibling seed's job, then this job's query step fails with
    # connection-refused since its own start_server_noninteractive.sh had already exited on that
    # mkdir error under set -e). Port 5001 is still hardcoded inside the .sif's flaskform.py (no
    # --port flag exists) and is NOT made seed-specific here -- a collision is only possible if
    # two seed jobs land on the very same compute node, which is not guaranteed to be safe.
    # SELENZYME2_DIR includes RHEA_TRACK too, not just SEED: within one job, this loop runs
    # Scaffold then Stratified sequentially, reusing the same seed -- a shared per-seed-only dir
    # let the Stratified iteration's REMAKE_SELENZYME2 rm-rf/rebuild race the still-shutting-down
    # Scaffold server (confirmed: Stratified silently produced zero results for two seeds --
    # query_selenzyme.py's query_reaction() swallows any non-200 server response and just skips
    # writing that reaction's file, no exception raised -- so the job "completed" with an empty
    # output directory and no visible error in the log beyond repeated "Processing reaction X...").
    DATA_SOURCE_DIR="$(pwd)/$FILTERED_DB_DIR" REMAKE_SELENZYME2=1 SELENZYME2_DIR="selenzyme2_seed${SEED}_${RHEA_TRACK}" \
        bash methods/SelenzymeRF/SelenzymeRF_code/start_server_noninteractive.sh &
    SERVER_PID=$!

    sleep 300

    RHEA_OUT_DIR="results/Case2/results-splits/seed_runs/seed${SEED}/${RHEA_TRACK}"
    mkdir -p "$RHEA_OUT_DIR"

    /scratch/jarcagniriv/Envs/clean_ec/bin/python3 methods/SelenzymeRF/SelenzymeRF_scripts/query_selenzyme.py \
        --csv_file "$RHEA_TEST_TSV" \
        --sep $'\t' \
        --results_folder "results/SelenzymeRF/case2_rhea_${RHEA_TRACK}_seed${SEED}" \
        --server_url http://localhost:5001 \
        --reaction_name_column REACTION_ID \
        --reaction_smiles_column REACTION_SMILES

    /scratch/jarcagniriv/Envs/clean_ec/bin/python3 methods/SelenzymeRF/SelenzymeRF_scripts/get_results.py \
        --folder "results/SelenzymeRF/case2_rhea_${RHEA_TRACK}_seed${SEED}" \
        --ec_column "EC Number" \
        --output "$RHEA_OUT_DIR/SelenzymeRF.csv"

    kill "$SERVER_PID" 2>/dev/null || true
done
