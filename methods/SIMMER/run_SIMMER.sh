#!/bin/bash
###################################################################################
### Run SIMMER
###################################################################################

#SBATCH --job-name=simmer
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jarcagniriv@unav.es
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH -o /scratch/jarcagniriv/Case2/SIMMER/logs/bench_%j.out

echo ===================================
echo ===     Load the Packages       ===
echo ===================================
echo `date`
module load Anaconda3
PYTHON=/scratch/jarcagniriv/Envs/SIMMERenv/bin/python

conda activate /Envs/SIMMERenv

########################################### Case 1 ############################################

# generate SIMMER input file
python /methods/SIMMER/SIMMER_scripts/simmer_input.py \
  --input-path /data/KEGG/kegg_reactions_with_split_can.csv \
  --output-path /methods/SIMMER/SIMMER_scripts/input/kegg_input_simmer.csv \
  --reaction-id-col Reaction ID \
  --reaction-sep "<=>" \
  --sp_col Equation \
  --smiles-col Isomeric_SMILES \
  --include-ec False

# query SIMMER 
python3 /SIMMER_code/SIMMER/SIMMER.py -i /SIMMER_code/SIMMER/SIMMER_files -o /methods/SIMMER/SIMMER_scripts/output/ResultsKEGG -q /methods/SIMMER/SIMMER_scripts/input/kegg_input_simmer.csv

python /methods/SIMMER/SIMMER_scripts/ec_predictions.py --input_dir /methods/SIMMER/SIMMER_scripts/output/ResultsKEGG --output_file results/Case1/SIMMER.csv
########################################### Case 2 ############################################

# generate SIMMER input file
python /methods/SIMMER/SIMMER_scripts/simmer_input.py \
  --input-path /data/MetaNetX/test_reactions.tsv \
  --output-path /methods/SIMMER/SIMMER_scripts/input/metanetx_input_simmer.csv \
  --reaction-id-col reaction_id \
  --sp_col substrates_products \
  --smiles-col reaction_smiles \
  --include-ec False

# generate SIMMER DB file
python /methods/SIMMER/SIMMER_scripts/simmer_input.py \
  --input-path /data/MetaNetX/train_reactions.tsv \
  --output-path /SIMMER_code/SIMMER/SIMMER_files_metanetx/chem_data/metanetx_reactions.csv \
  --reaction-id-col reaction_id \
  --sp_col substrates_products \
  --smiles-col reaction_smiles \
  --include-ec True \
  --ec-col ec

# generate SIMMER DB
python /methods/SIMMER/SIMMER_scripts/create_SIMMER_db.py
python /methods/SIMMER/SIMMER_scripts/ec_permutations.py

# query SIMMER
python3 /SIMMER_code/SIMMER/SIMMER2.py -i /SIMMER_code/SIMMER/SIMMER_files_metanetx -o /methods/SIMMER/SIMMER_scripts/output/ResultsMetaNetX -q /methods/SIMMER/SIMMER_scripts/input/metanetx_input_simmer.csv

python /methods/SIMMER/SIMMER_scripts/ec_predictions.py --input_dir /methods/SIMMER/SIMMER_scripts/output/ResultsMetaNetX --output_file results/Case2/SIMMER.csv

########################################### Case Study ############################################

# generate SIMMER input file
python /methods/SIMMER/SIMMER_scripts/simmer_input.py \
  --input-path /data/Drugs/drug_smiles_updated.csv \
  --output-path /methods/SIMMER/SIMMER_scripts/input/simmer_drugs_input.csv \
  --reaction-id-col drug \
  --sp_col right_comp \
  --smiles-col reaction_smiles \
  --include-ec False

# query SIMMER
python3 /SIMMER_code/SIMMER/SIMMER.py -i /SIMMER_code/SIMMER/SIMMER_files -o /methods/SIMMER/SIMMER_scripts/output/ResultsDrugs -q /methods/SIMMER/SIMMER_scripts/input/simmer_drugs_input.csv

python /methods/SIMMER/SIMMER_scripts/ec_predictions.py --input_dir /methods/SIMMER/SIMMER_scripts/output/ResultsDrugs --output_file results/CaseStudy/SIMMER.csv

