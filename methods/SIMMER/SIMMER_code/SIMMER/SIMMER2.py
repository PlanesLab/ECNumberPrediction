# SIMMER.py
# Annamarie Bustion
# 2022_02
#
# Example query:
#   $ python3 SIMMER.py -i <path_to>/SIMMER_files -o <path_to_output_dir>
######################################

#import mkl
import pickle as pkl
import pandas as pd
import numpy as np
from numpy import linalg as LA
import scipy.stats as stats
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import sys
import time
import argparse
import os
import shutil
import glob as glob

# Create reaction fingerprint for a query reaction.
def run_rxn(row, df, output_dir):
    sms_l = df.iloc[row, 3].split('.')
    sms_r = df.iloc[row, 4].split('.')
    
    smas_l = []
    for sm in sms_l:
        m = Chem.MolFromSmiles(sm.replace('R', '*'))
        sma = Chem.MolToSmarts(m, isomericSmiles=True)
        smas_l.append(sma)
    smas_r = []
    for sm in sms_r:
        m = Chem.MolFromSmiles(sm.replace('R', '*'))
        sma = Chem.MolToSmarts(m, isomericSmiles=True)
        smas_r.append(sma)
        
    left = '.'.join(map(str, smas_l))
    right = '.'.join(map(str, smas_r))
    
    rxn = rdChemReactions.ReactionFromSmarts(left + '>>' + right)
    
    # Reaction drawing (not needed for EC prediction)
    # drawer = rdMolDraw2D.MolDraw2DCairo(800, 200)
    # drawer.DrawReaction(rxn)
    # drawer.FinishDrawing()
    # drawer.WriteDrawingText(output_dir + '/' + df.iloc[row, 0] + '_input_reaction.png')
    
    return rxn

# Create reaction fingerprints for each query reaction.
def fp_queries(dms_df, fps, output_dir):
    print("\nDescribing input reactions...")
    t0 = time.time()
    for i in range(len(dms_df)):
        try:
            rxn = run_rxn(i, dms_df, output_dir)
            fp = Chem.rdChemReactions.CreateDifferenceFingerprintForReaction(rxn)
            fps.append(fp)
        except Exception as e:
            print("Improperly formatted SMILES; please check your input.", e)
            sys.exit()
    t1 = time.time()
    print("\nFinished creating fingerprints of queries in " + "{:.2f}".format(t1-t0) + " seconds")
    return fps

# Compute similarity matrix by adding query fingerprints to the precomputed Tanimoto matrix.
def add_queries_to_tanimoto(fps, tan_ar):
    t0 = time.time()
    metacyc_fps = fps[:34046]
    dms_fps = fps[34046:]
    tot_fps = metacyc_fps + dms_fps    
    tan_ar_dm = []
    for i in range(len(dms_fps)):
        tan = DataStructs.BulkTanimotoSimilarity(dms_fps[i], tot_fps)
        tan_ar_dm.append(tan)
    tot_array = np.append(np.array(tan_ar), np.array(tan_ar_dm)[:,:34046], axis=0)
    X_mol = np.append(tot_array, np.array(tan_ar_dm).T, axis=1)  
    t1 = time.time()
    print("\nFinished computing and adding queries to tanimoto similarity matrix in " 
          + "{:.2f}".format(t1-t0) + " seconds")
    return X_mol

# Find the closest reactions based on Euclidean distance in the similarity matrix.
def find_closest_rxns(query_id, X_mol, id_to_index):
    dists = {}
    vec = id_to_index[query_id]
    for i in range(len(X_mol)):
        dist = np.linalg.norm(X_mol[vec] - X_mol[i])
        dists[i] = dist
    s_dists = sorted(dists.items(), key=lambda x: x[1])
    results = []
    for ind, dist in s_dists:
        results.append([list(id_to_index.keys())[list(id_to_index.values()).index(ind)], str(dist)]) 
    [results.remove(result) for result in results if query_id in result]
    return [x for x in results if not x[0].startswith('DM')]

# Helper function for enrichment scoring.
def take_ES_walk(ec_cat, ec_df, level):
    running_tally = []
    i = 0
    for ec in ec_df[level]:
        if ec == ec_cat:
            i += 1
            running_tally.append(i)
        elif ec in ['NIL', 'DM']:
            running_tally.append(i)
        else:
            i -= 1
            running_tally.append(i)
    return running_tally

# Calculate enrichment p-values for EC predictions.
def find_EC_pval(ec_cats, ec_df, level, perm):
    sig_es_score_results = []
    for ec_cat in ec_cats:
        if ec_cat == 'NIL' or ec_cat == 'DM':
            continue
        else:
            es_tally = take_ES_walk(ec_cat, ec_df, level)
            es_score = max(es_tally)
            where = es_tally.index(es_score)
            es_score_norm = es_score / (ec_df[level].value_counts()[ec_cat])
            perm.append(es_score_norm)
            perm.sort(reverse=True)
            p_val = (perm.index(es_score_norm) + 1) / len(perm)
            perm.remove(es_score_norm)
            if p_val < 0.05 and es_score_norm > 0:
                sig_es_score_results.append([ec_cat, es_score_norm, p_val, where])
    return pd.DataFrame(sig_es_score_results, columns=['EC', 'enrich_score', 'p_val', 'where'])

# Predict EC numbers based on the ranked reaction list.
def predict_all_ECs(ranked_list, id_to_index, perm_df, rxn_to_ec):
    ec_array = []
    for result in ranked_list:
        ec_array.append([rxn_to_ec.get(result[0], 'NIL'), result[1]])
    ec_df = pd.DataFrame(ec_array, columns=['EC4', 'euc'])
    for ec in range(1, 4):
        ec_df['EC' + str(ec)] = ['.'.join(e.split('.')[0:ec]) for e in ec_df['EC4']]
    ec_df = ec_df[['euc', 'EC1', 'EC2', 'EC3', 'EC4']]
    
    # First level
    sig_results1 = find_EC_pval(ec_df['EC1'].unique(), ec_df, 'EC1', perm_df['EC1'].to_list()).sort_values(by='where')
    message = "a significant EC class prediction"
    
    if sig_results1.shape[0] > 0:
        ec_sub_df = ec_df.loc[ec_df['EC1'].isin(sig_results1['EC'])]
        sig_results2 = find_EC_pval(ec_sub_df['EC2'].unique(), ec_df, 'EC2', perm_df['EC2'].to_list()).sort_values(by='where')
        sig_results = pd.concat([sig_results1, sig_results2])
        
        if sig_results2.shape[0] > 0:
            ec_sub_df = ec_df.loc[ec_df['EC2'].isin(sig_results2['EC'])]
            sig_results3 = find_EC_pval(ec_sub_df['EC3'].unique(), ec_df, 'EC3', perm_df['EC3'].to_list()).sort_values(by='where')
            sig_results = pd.concat([sig_results, sig_results3])
            
            # Fourth level (extra code; left commented)
            # if sig_results3.shape[0] > 0:
            #     ec_sub_df = ec_df.loc[ec_df['EC3'].isin(sig_results3['EC'])]
            #     sig_results4 = find_EC_pval(ec_sub_df['EC4'].unique(), ec_df, 'EC4', perm_df['EC4'].to_list()).sort_values(by='where')
            #     sig_results = pd.concat([sig_results, sig_results4])
            # else:
            #     message = message + ", but no significant sub-sub-class"
        else:
            message = message + ", but no significant sub-class"
    else:
        message = "no significant EC prediction"
        sig_results = sig_results1
    
    return sig_results.reset_index(drop=True), message

# Return query results with only EC prediction functionality.
def return_query_results(DM, X_mol, id_to_index, rxn_to_ec, final, input_dir, output_dir,
                         perm_df, prot_dict, n_rxns):

    out_tsv = os.path.join(output_dir, f"{DM}_EC_predictions.tsv")
    if os.path.exists(out_tsv):
        print(f"[SKIP] {DM}: output already exists at {out_tsv}")
        return  
    
    t0 = time.time()
    ranked_list = find_closest_rxns(DM, X_mol, id_to_index)
    
    # EC predictions
    ecdf, message = predict_all_ECs(ranked_list, id_to_index, perm_df, rxn_to_ec)
    t1 = time.time()
    print("\nFor " + DM + ":\nFinished predicting EC codes in " + "{:.2f}".format(t1-t0) + " seconds,")
    print("Message: " + message)
    ecdf.to_csv(os.path.join(output_dir, DM + '_EC_predictions.tsv'), sep='\t', index=False)
    # Print predicted EC numbers using the 'EC' column from sig_results.
    print("\nPredicted EC numbers for query", DM, ":\n", ecdf['EC'].tolist())
    
    # Extra code (e.g., homolog predictions, distance ranking) is commented out.
    """
    # Homolog predictions code (commented)
    tsv_files = []
    fasta_files = []
    for i in range(n_rxns):
        num_lines = 0
        closest_rxn, euc_dist = ranked_list[i]
        try:
            files = prot_dict[closest_rxn.split('_')[0]]
        except KeyError:
            continue
        for file in files:
            f_tsv = glob.glob(input_dir + "/prot_data/" + file + '*tsv')
            num_lines = num_lines + sum(1 for line in open(f_tsv[0]))
            tsv_files.append(f_tsv[0])
            f_fasta = glob.glob(input_dir + "/prot_data/" + file + '*fasta')
            fasta_files.append(f_fasta[0])
        if num_lines < 1:
            print("Closest reaction is not associated with homologs. Try increasing n_rxns argument.")
    
    # Process TSV and FASTA files...
    """
    # End of extra code.

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_dir',
                        help='Input directory')
    parser.add_argument('-o', action='store', dest='output_dir',
                        help='Output directory')
    parser.add_argument('-q', action='store', dest='query',
                        help='Query tsv file location. Leave empty to input directly')
    parser.add_argument('-n', action='store', dest='n_rxns', default=1, type=int,
                        help='Number of reactions to query (default: 1)')
    # Additional arguments (threads, fingerprint type, PCA) commented out.
    return parser.parse_args()

def main():
    arguments = get_arguments()
    # mkl.set_num_threads(arguments.num_threads)
    input_dir = arguments.input_dir
    output_dir = arguments.output_dir
    query = arguments.query
    n_rxns = arguments.n_rxns

    # Load MetaCyc data.
    fps = pkl.load(open(os.path.join(input_dir, "chem_data/MC_rxn_fps.p"), 'rb'))
    tan_ar = pd.read_csv(os.path.join(input_dir, "chem_data/MC_rxn_tanimoto_matrix.csv"))
    id_to_index = pkl.load(open(os.path.join(input_dir, "chem_data/MC_rxn_to_matrix_index_dict.p"), 'rb'))
    rxn_to_ec = pkl.load(open(os.path.join(input_dir, "chem_data/MC_rxn_ec_dict.p"), 'rb'))
    final = pd.read_csv(os.path.join(input_dir, "chem_data/metanetx_reactions.csv"))
    perm_df = pd.read_csv(os.path.join(input_dir, "chem_data/ec_perm.csv"))
    prot_dict = pkl.load(open(os.path.join(input_dir, "prot_data/prot_dict.p"), 'rb'))
    print('\nPrecomputed data loaded')
    
    # Load and run query.
    if query is not None:
        dms_df = pd.read_csv(query, sep='\t')
        if dms_df.shape[1] != 5:
            dms_df = pd.read_csv(query, sep=',')
            if dms_df.shape[1] != 5:
                print("Improperly formatted input tsv. Make sure the file is tab or comma delimited with 5 columns.")
                sys.exit()
    else:
        reaction = input("\nEnter a reaction identifier (e.g. DM1): ")
        left_comp = input("\nEnter substrate name (e.g. gemcitabine): ")
        right_comp = input("\nEnter product name (e.g. 2′, 2′-difluorodeoxyuridine): ")
        left_smiles = input("\nEnter the SMILES of substrate(s) and cofactors separated by '.' : ") 
        right_smiles = input("\nEnter the SMILES of product(s) separated by '.' : ")
        if left_smiles or right_smiles:
            dms_df = pd.DataFrame([[reaction, left_comp, right_comp, left_smiles, right_smiles]],
                                  columns=['reaction', 'left_comp', 'right_comp', 'left_smiles', 'right_smiles'])
        else:
            print('No queries provided')
            sys.exit()
        
    X_mol = add_queries_to_tanimoto(fp_queries(dms_df, fps, output_dir), tan_ar)
    
    # Update indices for query reactions (assuming precomputed data has 8914 entries).
    for i in range(len(dms_df)):
        id_to_index[dms_df.iloc[i, 0]] = i + 34046
        rxn_to_ec[dms_df.iloc[i, 0]] = 'DM'
    
    # For each query reaction, predict EC numbers and output the list.
    for DM in dms_df['reaction'].values:
        return_query_results(DM, X_mol, id_to_index, rxn_to_ec, final, input_dir,
                             output_dir, perm_df, prot_dict, n_rxns)

if __name__ == '__main__':
    main()


