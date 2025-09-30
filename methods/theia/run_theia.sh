#!/bin/bash

###################################################################################
# Author: Josefina Arcagni
# Date: 9/9/2025
# Description: Run theia for various cases.
###################################################################################


#SBATCH --job-name=theia_GPU
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jarcagniriv@unav.es
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --gres=gpu:1      
#SBATCH -o /scratch/jarcagniriv/ECNumberPredictionReview/Theia/logs/bench_%j.out

module load Anaconda3  

python -c 'import sys; print(sys.version_info[:])'

# Optional: Display the PATH for debugging purposes
echo $PATH

################################################### CASE 1 ###########################################################################
conda activate theia

python /methods/theia/theia_scripts/query_theia.py --query_file /data/KEGG/8:2KEGGTest_canonicalized.txt 
    / --reaction_ids_file /data/KEGG/kegg_reactions_with_split_can.csv --reaction_id_column Reaction ID --output_file /results/Case1/theia.csv

################################################### CASE 2 ###########################################################################

conda activate theia

#create training DB and train the theia model
python /methods/theia/theia_code/theia/generateDB/encode_split_data.py --path /data/MetaNetX/train_metanetx.csv  --output-path /methods/theia/theia_code/theia/DB
./train_all.sh

#run theia
python /methods/theia/theia_scripts/query_theia.py --query_file /data/MetaNetX/queries.txt
    / --reaction_ids_file /data/MetaNetX/test_reactions.csv --reaction_id_column reaction_id --output_file /results/Case2/theia.csv --model DB.123

################################################### CASE STUDY ###########################################################################

conda activate theia

python /methods/theia/theia_scripts/query_theia.py --query_file /data/Drugs/reaction_smiles_can.txt
    / --reaction_ids_file /data/KEGG/kegg_reactions_with_split_can.csv --reaction_id_column drug --output_file /results/CaseStudy/theia.csv

