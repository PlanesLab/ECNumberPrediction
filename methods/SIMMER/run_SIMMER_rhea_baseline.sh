#!/bin/bash
###################################################################################
### Run SIMMER on the Rhea Stratified split (Case 2), cofactors LEFT IN.
### Baseline counterpart to run_SIMMER_rhea_nocofactor.sh, for direct comparison.
###################################################################################

#SBATCH --job-name=simmer_rhea_base
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jarcagniriv@unav.es
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH -o results/logs/simmer_rhea_base_%j.out

set -euo pipefail

echo "==================================="
echo "===     Load the Packages       ==="
echo "==================================="
echo "$(date)"
module load Miniforge3/24.11.3-2
eval "$(conda shell.bash hook)"
conda activate SIMMERenv
python3 -c "import rdkit" || { echo "FATAL: SIMMERenv did not activate correctly"; exit 1; }

SCRIPTS=methods/SIMMER/SIMMER_scripts
DB_DIR=methods/SIMMER/SIMMER_code/SIMMER/SIMMER_files_rhea_stratified_baseline
OUT_DIR=$SCRIPTS/output/ResultsRheaStratifiedBaseline
RESULTS_DIR=results/Case2/results-splits/SIMMER_baseline

mkdir -p "$OUT_DIR" "$RESULTS_DIR"

echo "==================================="
echo "===  1. Build reference DB      ==="
echo "==================================="
python3 "$SCRIPTS/build_rhea_nocofactor_db.py" \
    --input data/Splits-Rhea/Stratified/train.tsv \
    --output-dir "$DB_DIR"

echo "==================================="
echo "===  2. Build query file        ==="
echo "==================================="
python3 "$SCRIPTS/build_rhea_nocofactor_query.py" \
    --input data/Splits-Rhea/Stratified/test.tsv \
    --output-query "$SCRIPTS/input/rhea_stratified_baseline_query.csv" \
    --output-truth "$RESULTS_DIR/Stratified_baseline_ground_truth.csv"

echo "==================================="
echo "===  3. Run SIMMER2 (generic)   ==="
echo "==================================="
python3 methods/SIMMER/SIMMER_code/SIMMER/SIMMER2_generic.py \
    -i "$DB_DIR" \
    -o "$OUT_DIR" \
    -q "$SCRIPTS/input/rhea_stratified_baseline_query.csv"

echo "==================================="
echo "===  4. Aggregate predictions   ==="
echo "==================================="
python3 "$SCRIPTS/ec_predictions.py" \
    --input_dir "$OUT_DIR" \
    --output_file "$RESULTS_DIR/Stratified_SIMMER_baseline_predictions.csv"

echo "==================================="
echo "===  5. Evaluate (incl. class1) ==="
echo "==================================="
python3 "$SCRIPTS/evaluate_rhea_nocofactor.py" \
    --predictions "$RESULTS_DIR/Stratified_SIMMER_baseline_predictions.csv" \
    --ground-truth "$RESULTS_DIR/Stratified_baseline_ground_truth.csv" \
    --output "$RESULTS_DIR/Stratified_SIMMER_baseline_evaluation.csv"

echo "$(date)"
echo "Done."
