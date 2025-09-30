import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdChemReactions
import pickle
import os

# Define paths
input_path = "/SIMMER_code/SIMMER/SIMMER_files_metanetx/chem_data/metanetx_reactions.csv"
output_dir = "/SIMMER_code/SIMMER/SIMMER_files_metanetx/chem_data"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load the input database
try:
    input_db = pd.read_csv(input_path)
except Exception as e:
    raise FileNotFoundError(f"Could not load the input file: {e}")

# Filter out rows with empty EC numbers in the "EC_number" column.
def is_nonempty_ec(ec):
    ec = str(ec).strip()
    return ec != '' and ec.lower() != 'nan'

input_db = input_db[input_db['EC_number'].apply(is_nonempty_ec)]

# 1. Generate reaction fingerprints
def create_fingerprints(df):
    fps = []
    for i, row in df.iterrows():
        try:
            # Construct reaction SMARTS using left and right SMILES
            rxn = rdChemReactions.ReactionFromSmarts(f"{row['left_smiles']}>>{row['right_smiles']}")
            fp = rdChemReactions.CreateDifferenceFingerprintForReaction(rxn)
            fps.append(fp)
        except Exception as e:
            print(f"Error processing row {i}: {e}")
    return fps

reaction_fps = create_fingerprints(input_db)

# Save reaction fingerprints
with open(os.path.join(output_dir, "MC_rxn_fps.p"), "wb") as f:
    pickle.dump(reaction_fps, f)

# 2. Create a dictionary mapping reactions to EC numbers
def create_ec_dict(df):
    ec_dict = {}
    for _, row in df.iterrows():
        reaction = row["reaction"]
        ec_numbers = str(row["EC_number"]).split(',')
        # Use the first EC number if available
        ec_dict[reaction] = ec_numbers[0].strip() if ec_numbers else ""
    return ec_dict

reaction_ec_dict = create_ec_dict(input_db)

with open(os.path.join(output_dir, "MC_rxn_ec_dict.p"), "wb") as f:
    pickle.dump(reaction_ec_dict, f)

# 3. Create a Tanimoto similarity matrix as a normal pandas DataFrame
def create_tanimoto_matrix(fps):
    n = len(fps)
    matrix = np.zeros((n, n), dtype=np.float64)
    for i, fp in enumerate(fps):
        sims = DataStructs.BulkTanimotoSimilarity(fp, fps)
        matrix[i, :] = sims
        if i % 1000 == 0:
            print(f"Processed row {i} of {n}")
    return matrix

tanimoto_array = create_tanimoto_matrix(reaction_fps)
# Convert to pandas DataFrame
tanimoto_df = pd.DataFrame(tanimoto_array)

# Save the Tanimoto similarity matrix to CSV
tanimoto_matrix_file = os.path.join(output_dir, "MC_rxn_tanimoto_matrix.csv")
tanimoto_df.to_csv(tanimoto_matrix_file, index=False)

# 4. Create a mapping of reactions to matrix indices
def create_index_mapping(df):
    return {reaction: idx for idx, reaction in enumerate(df["reaction"])}

reaction_to_index = create_index_mapping(input_db)

with open(os.path.join(output_dir, "MC_rxn_to_matrix_index_dict.p"), "wb") as f:
    pickle.dump(reaction_to_index, f)

# 5. Define functions for EC enrichment and p-value calculation
def take_ES_walk(ec_cat, ec_df, level):
    running_tally = []
    i = 0
    for ec in ec_df[level]:
        if ec == ec_cat:
            i += 1
        elif ec in ['NIL', 'DM']:
            pass  # No change for NIL or DM
        else:
            i -= 1
        running_tally.append(i)
    return running_tally

def find_EC_pval(ec_cats, ec_df, level, num_permutations=1000):
    all_results = []
    for ec_cat in ec_cats:
        if ec_cat in ['NIL', 'DM']:
            continue  # Skip NIL and DM categories
        # Calculate enrichment score for the actual data
        es_tally = take_ES_walk(ec_cat, ec_df, level)
        es_score = max(es_tally)
        where = es_tally.index(es_score)
        es_score_norm = es_score / ec_df[level].value_counts().get(ec_cat, 1)
        # Generate permutation scores
        permuted_scores = []
        for _ in range(num_permutations):
            permuted_labels = ec_df[level].sample(frac=1, replace=False).tolist()
            perm_es_tally = take_ES_walk(ec_cat, pd.DataFrame({level: permuted_labels}), level)
            perm_es_score = max(perm_es_tally) / ec_df[level].value_counts().get(ec_cat, 1)
            permuted_scores.append(perm_es_score)
        permuted_scores = sorted(permuted_scores, reverse=True)
        p_val = sum(perm_score >= es_score_norm for perm_score in permuted_scores) / num_permutations
        all_results.append([level, ec_cat, es_score_norm, p_val, where])
    return all_results

def compute_ec_levels(ec_df, num_permutations=1000):
    final_results = pd.DataFrame()
    levels = ['EC1', 'EC2', 'EC3', 'EC4']
    level_names = [1, 2, 3, 4]
    for level_name, level_col in zip(level_names, levels):
        ec_df[level_col] = ec_df['EC'].apply(lambda x: '.'.join(x.split('.')[:level_name]) if pd.notna(x) else x)
        ec_cats = ec_df[level_col].unique()
        ec_pvals_df = find_EC_pval(ec_cats, ec_df, level_col, num_permutations)
        ec_pvals_df = pd.DataFrame(ec_pvals_df, columns=['Level', 'EC', 'enrich_score', 'p_val', 'where'])
        ec_pvals_df['Level'] = level_col
        final_results = pd.concat([final_results, ec_pvals_df], ignore_index=True)
    return final_results

# 6. Calculate p-values and save results
# Build a DataFrame from the reaction-to-EC dictionary
ec_df = pd.DataFrame.from_dict(reaction_ec_dict, orient='index', columns=['EC']).reset_index().rename(columns={'index': 'Reaction'})
final_pval_results = compute_ec_levels(ec_df)

# Save the final results (uncomment the next line to actually write the CSV)
final_pval_file = os.path.join(output_dir, "final_ec_pvals.csv")
# final_pval_results.to_csv(final_pval_file, index=False)

print("Processing complete. Outputs saved in:")
print(f"1. {os.path.join(output_dir, 'MC_rxn_fps.p')} (Reaction Fingerprints)")
print(f"2. {os.path.join(output_dir, 'MC_rxn_ec_dict.p')} (Reaction-to-EC Dictionary)")
print(f"3. {tanimoto_matrix_file} (Tanimoto similarity matrix)")
print(f"4. {os.path.join(output_dir, 'MC_rxn_to_matrix_index_dict.p')} (Reaction-to-Index Mapping)")
print(f"5. {final_pval_file} (Final EC P-Values)")
