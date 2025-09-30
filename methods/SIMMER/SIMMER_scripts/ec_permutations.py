import pandas as pd

def take_ES_walk(ec_cat, ec_df, level):
    """
    Walk over the EC values for the given level.
    For each row, add +1 if the value matches ec_cat,
    subtract 1 if it does not (except for 'NIL' or 'DM', which are ignored).
    Returns a list representing the running tally.
    """
    running_tally = []
    tally = 0
    for ec in ec_df[level]:
        if ec == ec_cat:
            tally += 1
            running_tally.append(tally)
        elif ec in ['NIL', 'DM']:
            running_tally.append(tally)
        else:
            tally -= 1
            running_tally.append(tally)
    return running_tally

def compute_perm_score(ec_df, level, ec_cat):
    """
    Shuffles the values in the specified level column,
    computes the enrichment walk for the given ec_cat,
    and returns the normalized maximum score.
    """
    # Shuffle the column values
    permuted_labels = ec_df[level].sample(frac=1, replace=False).tolist()
    # Create a temporary DataFrame using the permuted labels
    temp_df = pd.DataFrame({ level: permuted_labels })
    # Compute the enrichment walk for the representative EC
    es_tally = take_ES_walk(ec_cat, temp_df, level)
    es_score = max(es_tally)
    # Normalize by the count of ec_cat in the original data
    denom = ec_df[level].value_counts().get(ec_cat, 1)
    return es_score / denom

def main(num_permutations=1000, output_csv='permutation_scores.csv'):
    # Read the CSV file.
    file_path = '/SIMMER_code/SIMMER/SIMMER_files_metanetx/chem_data/metanetx_reactions.csv'
    df = pd.read_csv(file_path)
    
    # Ensure the EC_number column is treated as a string.
    df['EC_number'] = df['EC_number'].astype(str)
    
    # If a row has more than one EC number (separated by '|'),
    # split them into separate rows.
    df = df.assign(EC_number=df['EC_number'].str.split('|')).explode('EC_number')
    # Remove any extra whitespace around the EC numbers.
    df['EC_number'] = df['EC_number'].str.strip()
    
    # Create hierarchical EC columns.
    df['EC1'] = df['EC_number'].apply(lambda x: x.split('.')[0] if isinstance(x, str) else None)
    df['EC2'] = df['EC_number'].apply(lambda x: '.'.join(x.split('.')[0:2]) if isinstance(x, str) else None)
    df['EC3'] = df['EC_number'].apply(lambda x: '.'.join(x.split('.')[0:3]) if isinstance(x, str) else None)
    df['EC4'] = df['EC_number']  # Full EC number
    
    # Create a DataFrame with just the EC columns (drop rows with missing data)
    ec_df = df[['EC1', 'EC2', 'EC3', 'EC4']].dropna()
    
    # For each level, remove special values and choose a representative EC category.
    valid_EC1 = ec_df[~ec_df['EC1'].isin(['NIL', 'DM'])]
    valid_EC2 = ec_df[~ec_df['EC2'].isin(['NIL', 'DM'])]
    valid_EC3 = ec_df[~ec_df['EC3'].isin(['NIL', 'DM'])]
    valid_EC4 = ec_df[~ec_df['EC4'].isin(['NIL', 'DM'])]
    
    rep_EC1 = valid_EC1['EC1'].mode()[0] if not valid_EC1.empty else None
    rep_EC2 = valid_EC2['EC2'].mode()[0] if not valid_EC2.empty else None
    rep_EC3 = valid_EC3['EC3'].mode()[0] if not valid_EC3.empty else None
    rep_EC4 = valid_EC4['EC4'].mode()[0] if not valid_EC4.empty else None
    
    print("Representative EC categories selected:")
    print("EC1:", rep_EC1)
    print("EC2:", rep_EC2)
    print("EC3:", rep_EC3)
    print("EC4:", rep_EC4)
    
    # Compute permutation scores for each level for each iteration.
    perm_scores_list = []
    for _ in range(num_permutations):
        score_EC1 = compute_perm_score(ec_df, 'EC1', rep_EC1) if rep_EC1 is not None else None
        score_EC2 = compute_perm_score(ec_df, 'EC2', rep_EC2) if rep_EC2 is not None else None
        score_EC3 = compute_perm_score(ec_df, 'EC3', rep_EC3) if rep_EC3 is not None else None
        score_EC4 = compute_perm_score(ec_df, 'EC4', rep_EC4) if rep_EC4 is not None else None
        
        perm_scores_list.append([score_EC1, score_EC2, score_EC3, score_EC4])
    
    # Create a DataFrame for the permutation results.
    perm_df = pd.DataFrame(perm_scores_list, columns=['EC1', 'EC2', 'EC3', 'EC4'])
    
    # Output the permutation scores to a CSV file.
    perm_df.to_csv(output_csv, index=False)
    print(f"Permutation scores saved to {output_csv}")
    print(perm_df.head())

if __name__ == '__main__':
    main(num_permutations=1000, output_csv='/SIMMER_code/SIMMER/SIMMER_files_metanetx/chem_data/ec_perm.csv')
