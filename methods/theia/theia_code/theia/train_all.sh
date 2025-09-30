SCRIPT_PATH='/theia/scripts'
DATA_PATH='/theia/DB'
MODEL_PATH='/theia/models'

python ${SCRIPT_PATH}/train.py ${DATA_PATH}/DB-0-ec123-train.csv ${DATA_PATH}/DB-0-ec123-valid.csv ${DATA_PATH}/DB-0-ec123-test.csv ${MODEL_PATH}/DB-0-ec123
