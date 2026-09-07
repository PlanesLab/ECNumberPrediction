#!/bin/bash
# Usage: train_all.sh <data_path> <model_path> <prefix> [seed]
# e.g.:  train_all.sh /scratch/.../DB /scratch/.../models DB-0-ec123 0
set -euo pipefail

DATA_PATH="$1"
MODEL_PATH="$2"
PREFIX="$3"
SEED="${4:-42}"
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")/scripts" && pwd)"

mkdir -p "$MODEL_PATH"

# NOTE: train.py's @click.argument order is (train, test, valid, output) --
# not (train, valid, test, output) despite the function signature reading
# that way -- click binds CLI positionals to argument names in decorator
# order, so passing valid/test swapped here would silently feed the model
# selection step the real test set. Keep this order matching the decorators.
python "${SCRIPT_PATH}/train.py" \
    "${DATA_PATH}/${PREFIX}-train.csv" \
    "${DATA_PATH}/${PREFIX}-test.csv" \
    "${DATA_PATH}/${PREFIX}-valid.csv" \
    "${MODEL_PATH}/${PREFIX}" \
    --seed "${SEED}"
