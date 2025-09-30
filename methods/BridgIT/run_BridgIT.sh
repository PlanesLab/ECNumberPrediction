# -----------------------------------------------------------------------------
# Author: jarcagniriv
# Script: run_BridgIT.sh
# Description: Prepares input for the BridgIT web server and extracts EC numbers from the output.
# -----------------------------------------------------------------------------

# Case 1 input
python /methods/BridgIT/bridgIT_scripts/bridgit_input.py --input_file /data/KEGG/kegg_reactions_current_test.csv 
/ --output_file /methods/BridgIT/input/reduced_inputs_KEGG/systemfile.txt --molfile_folder /methods/BridgIT/input/reduced_inputs_KEGG/molfiles
/ --equation_column Equation

# The BridgIT server is slow in the case of a large number of queries. Therefore it is recommended to split the input file into smaller chunks.
# in that case use the /methods/BridgIT/bridgIT_scripts/bridgit_reduced_input.py script to do so.

# Case Study Input
python /methods/BridgIT/bridgIT_scripts/bridgit_input_drugs.py

# Input the zip files in: https://lcsb-databases.epfl.ch/Bridgit (User required for use)

# After receiving the results from the BridgIT server, extract the EC numbers using the following script:

# Case 1 

python /methods/BridgIT/bridgIT_scripts/get_results.py --bridgit_dir /methods/BridgIT/input/reduced_inputs_KEGG --input_csv /data/KEGG/kegg_reactions_current_test.csv
/ --output /results/Case1/BridgIT.csv --csv_reaction_col Reaction ID --csv_ec_col EC Number
# Case Study

python /methods/BridgIT/bridgIT_scripts/get_results.py --bridgit_dir /methods/BridgIT/input/input_casestudy --input_csv /data/Drugs/drug_smiles_updated.csv
/ --output /results/CaseStudy/BridgIT.csv --csv_reaction_col drug --csv_ec_col ec


