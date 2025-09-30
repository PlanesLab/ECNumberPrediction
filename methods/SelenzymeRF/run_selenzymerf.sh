#!/bin/bash
###################################################################################
# Author: Josefina Arcagni
# Date: 9/9/2025
# Description: Run SelenzymeRF for various cases.
###################################################################################

###################################################################################

echo "Running SelenzymeRF for Case 1"
bash /SelenzymeRF_code/start_server.sh &

# Wait for the server to start (5 minutes)
sleep 300

# Run the query script
python3 /SelenzymeRF_scripts/query_selenzyme.py --csv_file /data/KEGG/kegg_reactions_current_test.csv --results_folder /results/SelenzymeRF_case1 --server_url http://localhost:5000 --reaction_name_column Reaction ID --reaction_smiles_column Isomeric_SMILES

# Run the results script
python3 /SelenzymeRF_scripts/get_results.py --folder /results/SelenzymeRF_case1 --reaction_csv /data/KEGG/kegg_reactions_current_test.csv --reaction_column Reaction ID --ec_column EC Number --output /results/Case1/SelenzymeRF.csv

######################################################################################################################################################################################################################################################

echo "Running SelenzymeRF for Case Study"
bash /SelenzymeRF_code/start_server.sh &

# Wait for the server to start (5 minutes)
sleep 300

# Run the query script
python3 /SelenzymeRF_scripts/query_selenzyme.py --csv_file /data/Drugs/drug_smiles_updated.csv --results_folder /results/SelenzymeRF_casestudy --server_url http://localhost:5000 --reaction_name_column drug --reaction_smiles_column reaction_smiles

# Run the results script
python3 /SelenzymeRF_scripts/get_results.py --folder /results/SelenzymeRF_casestudy --reaction_csv /data/Drugs/drug_smiles_updated.csv --reaction_column drug --ec_column EC Number --output /results/CaseStudy/SelenzymeRF.csv

################################################################################################################################################################

echo "Running SelenzymeRF for Case 2"

# Modify database, replace new one with original one
python3 /methods/SelenzymeRF/SelenzymeRF_scripts/generate_selenzyme_db.py

bash /SelenzymeRF_code/start_server.sh &

# Wait for the server to start (5 minutes)
sleep 300

# Run the query script
python3 /SelenzymeRF_scripts/query_selenzyme.py --csv_file /data/Drugs/drug_smiles_updated.csv --results_folder /results/SelenzymeRF_casestudy --server_url http://localhost:5000 --reaction_name_column drug --reaction_smiles_column reaction_smiles

# Run the results script
python3 /SelenzymeRF_scripts/get_results.py --folder /results/SelenzymeRF_casestudy --reaction_csv /data/Drugs/drug_smiles_updated.csv --reaction_column drug --ec_column EC Number --output /results/Case2/SelenzymeRF.csv