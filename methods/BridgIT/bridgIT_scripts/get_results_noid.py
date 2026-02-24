import os
import zipfile
import pandas as pd
import argparse
import re

def extract_all_ec_predictions(file_path, ec_col):
    try:
        df = pd.read_csv(file_path, sep="\t")
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return ""
    if df.empty or ec_col not in df.columns:
        print(f"Missing '{ec_col}' in {file_path}. Columns found: {df.columns.tolist()}")
        return ""
    ec_entries = []
    seen_prefixes = set()
    for val in df[ec_col]:
        s = str(val).strip().strip(";")
        parts = s.split("/")
        if len(parts) >= 3:
            ec_str = parts[2].strip()
            if ec_str:
                for ec in [e.strip() for e in ec_str.split(",") if e.strip()]:
                    prefix = ".".join(ec.split(".")[:3])
                    if prefix not in seen_prefixes:
                        seen_prefixes.add(prefix)
                        ec_entries.append(ec)
    return ";".join(ec_entries)

def extract_reaction_id_from_filename(filename):
    """
    Extract reaction ID from filename like 'Tanimoto_Atom_reactions_42.txt'
    Returns '42' in this example
    """
    # Remove .txt extension
    name = filename.replace('.txt', '')
    
    # Try to find the last number in the filename
    match = re.search(r'_(\d+)$', name)
    if match:
        return match.group(1)
    
    # Fallback: return the whole filename without extension
    return name

def process_zip_and_create_df(zip_path, extract_dir, ec_col):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)

    data = []
    for root, dirs, files in os.walk(extract_dir):
        for file in files:
            if file.endswith('.txt'):
                file_path = os.path.join(root, file)
                ec_pred = extract_all_ec_predictions(file_path, ec_col)
                
                # Extract reaction ID from filename
                reaction_id = extract_reaction_id_from_filename(file)
                
                data.append({
                    'reaction': reaction_id,
                    'EC_number_predicted': ec_pred
                })
    return pd.DataFrame(data)

def process_bridgit_results(results_dir, ec_col):
    base = os.path.join(results_dir, "temp_extracted")
    all_dfs = []
    for fname in os.listdir(results_dir):
        if fname.lower().endswith(".zip"):
            zipf = os.path.join(results_dir, fname)
            ext = os.path.join(base, os.path.splitext(fname)[0])
            os.makedirs(ext, exist_ok=True)
            df = process_zip_and_create_df(zipf, ext, ec_col)
            all_dfs.append(df)
    if all_dfs:
        return pd.concat(all_dfs, ignore_index=True)
    return pd.DataFrame(columns=["reaction", "EC_number_predicted"])

def main():
    parser = argparse.ArgumentParser(description="Process BridgIT results extracting reaction IDs from filenames")
    parser.add_argument("-b", "--bridgit_dir", required=True, help="Directory containing BridgIT zip results")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file path")
    parser.add_argument("--bridgit_ec_col", default="reactionsA/ECA", help="Column in BridgIT txt files with EC predictions")
    args = parser.parse_args()

    bridgit_df = process_bridgit_results(args.bridgit_dir, args.bridgit_ec_col)
    
    # Remove duplicates and save
    final = bridgit_df[['reaction', 'EC_number_predicted']].drop_duplicates()
    final.to_csv(args.output, index=False)
    print(f"Saved {len(final)} results to {args.output}")

if __name__ == "__main__":
    main()