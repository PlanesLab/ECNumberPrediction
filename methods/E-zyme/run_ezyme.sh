# -----------------------------------------------------------------------------
# Author: jarcagniriv
# Script: run_Ezyme.sh
# Description: E-zyme 1 and 2 webscrapping and results for EC number prediction.
# -----------------------------------------------------------------------------

# Case 1 webscrapping (inputs are Substrate-Product pairs from KEGG)
python /methods/E-zyme/ezyme_scripts/ezyme_webscrapping.py \
    -i /scratch/jarcagniriv/ECNumberPredictionReview/E-zyme/input/SubstrateProductPairs82.csv \
    -o /methods/E-zyme/output/outputKEGG \
    --reaction_id_col "Reaction ID" \
    --reactant_col "Reactants" \
    --product_col "Products" \
    --delimiter ";"

# Case Study webscrapping (inputs are manually curated Substrate-Product pairs from drugs)
python /methods/E-zyme/ezyme_scripts/ezyme_webscrapping.py  \
    -i /scratch/jarcagniriv/CaseStudy/drugs/sp_pairs_drugs.csv \
    -o /methods/E-zyme/output/outputDrugs \
    --drug_col drug \
    --pair1_col Pair1 \
    --pair2_col Pair2 \
    --delimiter ";"

# After receiving the results from the E-zyme server, extract the EC numbers using the following script:º

python /methods/E-zyme/ezyme_scripts/get_ezyme_results.py --input_dir /methods/E-zyme/output/outputKEGG --output_file results/Case1/E-zyme.csv

python /methods/E-zyme/ezyme_scripts/get_ezyme_results.py --input_dir /methods/E-zyme/output/outputDrugs --output_file results/CaseStudy/E-zyme.csv