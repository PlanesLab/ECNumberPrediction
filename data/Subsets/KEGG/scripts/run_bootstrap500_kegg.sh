#!/bin/bash
# Orchestrates the 10x500 KEGG bootstrap for one method: for each seed 0-9,
# splits off that seed's uncached ("todo") reactions via kegg_cache.py, and if
# non-empty, submits that method's existing run_X.sh pointed at the todo subset
# (via the SEED="<n>_<Method>" string-suffix trick -- run_X.sh's Case 1 block
# already does `if [ -d "data/Subsets/KEGG/seed_pools/seed${SEED}" ]`, so any
# string works, not just a number).
#
# Does NOT run the merge step (that must happen after the SLURM job finishes) --
# run merge_bootstrap500_kegg.sh once jobs complete.
#
# Usage: run_bootstrap500_kegg.sh <Method> <run_script_path> [extra sbatch args...]
set -euo pipefail

METHOD="$1"
RUN_SCRIPT="$2"
shift 2
EXTRA_SBATCH_ARGS=("$@")

PYTHON=/scratch/jarcagniriv/Envs/clean_ec/bin/python3
SAMPLES_DIR="data/Subsets/KEGG/bootstrap500_samples"
LOG="/tmp/claude-6347/-scratch-jarcagniriv/bab2e5c2-ab50-4ce8-9eb7-6b0699796605/scratchpad/bootstrap500_submitted_${METHOD}.tsv"
echo -e "seed\tjob_id\ttodo_count" > "$LOG"

for SEED in 0 1 2 3 4 5 6 7 8 9; do
    SEED_TAG="${SEED}_kegg500${METHOD}"
    TODO_DIR="data/Subsets/KEGG/seed_pools/seed${SEED_TAG}"

    $PYTHON data/Subsets/KEGG/scripts/kegg_cache.py split \
        --method "$METHOD" \
        --seed_sample "$SAMPLES_DIR/seed${SEED}/kegg_reactions_current_test.csv" \
        --output_dir "$TODO_DIR" > /tmp/split_${METHOD}_${SEED}.log 2>&1

    TODO_COUNT=$(($(wc -l < "$TODO_DIR/kegg_reactions_current_test.csv") - 1))
    echo "Seed $SEED ($METHOD): $TODO_COUNT to query"

    if [ "$TODO_COUNT" -le 0 ]; then
        echo -e "${SEED}\tSKIPPED_EMPTY\t0" >> "$LOG"
        continue
    fi

    OUT=$(sbatch --job-name="${METHOD}_kegg500_seed${SEED}" "${EXTRA_SBATCH_ARGS[@]}" \
        --export=ALL,SEED=${SEED_TAG},SKIP_METANETX=1,SKIP_CASESTUDY=1 \
        "$RUN_SCRIPT" 2>&1)
    JOBID=$(echo "$OUT" | grep -oE '[0-9]+$')
    echo "$OUT"
    echo -e "${SEED}\t${JOBID}\t${TODO_COUNT}" >> "$LOG"
done
