import os
import pandas as pd
import argparse

def extract_subsubclass(ec_number):
    """Extract the subsubclass (first 3 levels) from an EC number.
    Example: '1.2.3.4' -> '1.2.3'
    """
    parts = str(ec_number).split('.')
    if len(parts) >= 3:
        return '.'.join(parts[:3])
    return str(ec_number)

def collapse_duplicates_in_rank(rank_group):
    """Collapse duplicate subsubclasses within a single rank.
    
    Input: "2.4.1.1|2.4.1.2|2.4.1.3|3.1.1.1"
    Output: "2.4.1|3.1.1"
    """
    if not rank_group:
        return ""
    
    predictions = rank_group.split('|')
    seen_subsubclasses = set()
    collapsed = []
    
    for pred in predictions:
        pred = pred.strip()
        if not pred or pred == 'nan' or pred == 'NaN':
            continue
        subsubclass = extract_subsubclass(pred)
        if subsubclass not in seen_subsubclasses:
            seen_subsubclasses.add(subsubclass)
            collapsed.append(subsubclass)
    
    return '|'.join(collapsed) if collapsed else ""

def collapse_predictions(prediction_string):
    """Collapse duplicate subsubclasses in each ranking level.
    
    Input: "2.4.1.1|2.4.1.2;3.1.1.1|3.1.1.2"
    Output: "2.4.1;3.1.1"
    """
    if not prediction_string:
        return ""
    
    ranks = prediction_string.split(';')
    collapsed_ranks = []
    
    for rank in ranks:
        collapsed = collapse_duplicates_in_rank(rank)
        if collapsed:
            collapsed_ranks.append(collapsed)
    
    return ';'.join(collapsed_ranks)

def interleave_rankings(predictions_list):
    """Interleave rankings from multiple pairs.
    
    Input: List of prediction strings, each with rankings separated by semicolons
           e.g., ["1.2.3.4;1.2.3.5", "2.3.4.5;2.3.4.6"]
    
    Output: Single string with interleaved rankings
           e.g., "1.2.3.4|2.3.4.5;1.2.3.5|2.3.4.6"
    """
    if not predictions_list:
        return ""
    
    # Split each prediction string into individual rankings
    all_rankings = [pred.split(';') for pred in predictions_list]
    
    # Find maximum ranking length
    max_length = max(len(rankings) for rankings in all_rankings)
    
    # Interleave rankings by position
    interleaved = []
    for i in range(max_length):
        rank_group = []
        for rankings in all_rankings:
            if i < len(rankings):
                rank_group.append(rankings[i])
        if rank_group:
            interleaved.append('|'.join(rank_group))
    
    return ';'.join(interleaved)

def process_folder(folder_path, many_pairs, debug=False):
    """Process a single reaction folder and extract E-zyme predictions.
    
    Args:
        folder_path: Path to the reaction folder
        many_pairs: Boolean indicating whether to look for pair subdirectories
        debug: If True, print debug information
    
    Returns:
        Tuple of (ezyme1_result, ezyme2_result)
    """
    ezyme1_predictions = []
    ezyme2_predictions = []
    
    if many_pairs:
        # Look for pair_* subdirectories
        items = os.listdir(folder_path)
        pair_folders = sorted([item for item in items if item.startswith('pair_') and 
                               os.path.isdir(os.path.join(folder_path, item))])
        
        if debug:
            print(f"\n[DEBUG] Processing {folder_path}")
            print(f"[DEBUG] Found {len(pair_folders)} pair folders: {pair_folders}")
        
        if pair_folders:
            for pair_folder in pair_folders:
                pair_path = os.path.join(folder_path, pair_folder)
                e1_preds, e2_preds = extract_predictions_from_path(pair_path, debug=debug)
                if debug:
                    print(f"[DEBUG] {pair_folder}: E-zyme1={e1_preds}, E-zyme2={e2_preds}")
                ezyme1_predictions.extend(e1_preds)
                ezyme2_predictions.extend(e2_preds)
        else:
            if debug:
                print(f"[DEBUG] No pair folders found, processing folder directly")
            e1_preds, e2_preds = extract_predictions_from_path(folder_path, debug=debug)
            ezyme1_predictions.extend(e1_preds)
            ezyme2_predictions.extend(e2_preds)
    else:
        e1_preds, e2_preds = extract_predictions_from_path(folder_path, debug=debug)
        ezyme1_predictions.extend(e1_preds)
        ezyme2_predictions.extend(e2_preds)
    
    if debug:
        print(f"[DEBUG] Before interleaving: E-zyme1={ezyme1_predictions}, E-zyme2={ezyme2_predictions}")
    
    # Interleave rankings from different pairs
    ezyme1_result = interleave_rankings(ezyme1_predictions)
    ezyme2_result = interleave_rankings(ezyme2_predictions)
    
    if debug:
        print(f"[DEBUG] After interleaving: E-zyme1={ezyme1_result}, E-zyme2={ezyme2_result}")
    
    # Collapse duplicate subsubclasses in each rank
    ezyme1_result = collapse_predictions(ezyme1_result)
    ezyme2_result = collapse_predictions(ezyme2_result)
    
    if debug:
        print(f"[DEBUG] After collapsing: E-zyme1={ezyme1_result}, E-zyme2={ezyme2_result}")
    
    return ezyme1_result, ezyme2_result

def extract_predictions_from_path(path, debug=False):
    """Extract E-zyme predictions from CSV files in a given path.
    
    For each file, collect all predictions (rows) and join them with semicolons.
    Returns a single string for each file found, with rankings separated by semicolons.
    
    Args:
        path: Directory path to search for E-zyme CSV files
        debug: If True, print debug information
    
    Returns:
        Tuple of (ezyme1_predictions, ezyme2_predictions) as lists of strings
    """
    ezyme1_predictions = []
    ezyme2_predictions = []
    
    try:
        files = os.listdir(path)
        if debug:
            print(f"[DEBUG]   Files in {path}: {files}")
    except Exception as e:
        print(f"Error listing directory {path}: {e}")
        return ezyme1_predictions, ezyme2_predictions
    
    for file in files:
        # Check for both naming conventions: EZYME1.csv and *_E-zyme1.csv
        if file == 'EZYME1.csv' or '_E-zyme1.csv' in file:
            ezyme1_file = os.path.join(path, file)
            try:
                ezyme1_df = pd.read_csv(ezyme1_file)
                if debug:
                    print(f"[DEBUG]   E-zyme1 file: {file}, shape: {ezyme1_df.shape}")
                    if not ezyme1_df.empty:
                        print(f"[DEBUG]   First few values: {ezyme1_df.iloc[:3, 0].tolist()}")
                if not ezyme1_df.empty and ezyme1_df.shape[0] > 0:
                    preds = ezyme1_df.iloc[:, 0].astype(str).str.strip().tolist()
                    preds = [p for p in preds if p and p != 'nan' and p != 'NaN']
                    if preds:
                        pred_string = ';'.join(preds)
                        ezyme1_predictions.append(pred_string)
                        if debug:
                            print(f"[DEBUG]   Added E-zyme1 predictions: {pred_string}")
            except Exception as e:
                print(f"Error reading E-zyme1 file {ezyme1_file}: {e}")
        
        # Check for both naming conventions: EZYME2.csv and *_E-zyme2.csv
        elif file == 'EZYME2.csv' or '_E-zyme2.csv' in file:
            ezyme2_file = os.path.join(path, file)
            try:
                ezyme2_df = pd.read_csv(ezyme2_file)
                if debug:
                    print(f"[DEBUG]   E-zyme2 file: {file}, shape: {ezyme2_df.shape}")
                    if not ezyme2_df.empty and ezyme2_df.shape[1] > 2:
                        print(f"[DEBUG]   First few values: {ezyme2_df.iloc[:3, 2].tolist()}")
                if not ezyme2_df.empty and ezyme2_df.shape[1] > 2 and ezyme2_df.shape[0] > 0:
                    preds = ezyme2_df.iloc[:, 2].astype(str).str.strip().tolist()
                    preds = [p for p in preds if p and p != 'nan' and p != 'NaN']
                    if preds:
                        pred_string = ';'.join(preds)
                        ezyme2_predictions.append(pred_string)
                        if debug:
                            print(f"[DEBUG]   Added E-zyme2 predictions: {pred_string}")
            except Exception as e:
                print(f"Error reading E-zyme2 file {ezyme2_file}: {e}")
    
    return ezyme1_predictions, ezyme2_predictions

def main():
    parser = argparse.ArgumentParser(description="Collect E-zyme predictions from multiple folders.")
    parser.add_argument("--input_dir", required=True, 
                        help="Path to the root directory containing reaction folders")
    parser.add_argument("--output_file", required=True, 
                        help="Path to save the final CSV file")
    parser.add_argument("--many_pairs", type=lambda x: x.lower() == 'true', 
                        default=False,
                        help="If true, look for pair_* subdirectories within each reaction folder (default: False)")
    parser.add_argument("--debug", action='store_true',
                        help="Enable debug output")
    args = parser.parse_args()
    
    root_dir = args.input_dir
    csv_data = []
    
    for idx, reaction_folder in enumerate(sorted(os.listdir(root_dir))):
        reaction_folder_path = os.path.join(root_dir, reaction_folder)
        
        if os.path.isdir(reaction_folder_path):
            debug = args.debug and idx < 3
            ezyme1_result, ezyme2_result = process_folder(reaction_folder_path, args.many_pairs, debug=debug)
            csv_data.append([reaction_folder, ezyme1_result, ezyme2_result])
    
    # Save final CSV
    final_df = pd.DataFrame(csv_data, columns=['Reaction Name', 'E-zyme1', 'E-zyme2'])
    final_df.to_csv(args.output_file, index=False)
    print(f"CSV file created successfully: {args.output_file}")
    print(f"Processed {len(csv_data)} reaction folders")

if __name__ == "__main__":
    main()