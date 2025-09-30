import os
import zipfile
import pandas as pd
import argparse

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

def process_zip_and_create_df(zip_path, extract_dir, ec_col):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)

    data = []
    for root, dirs, files in os.walk(extract_dir):
        for file in files:
            if file.endswith('.txt'):
                file_path = os.path.join(root, file)
                ec_pred = extract_all_ec_predictions(file_path, ec_col)
                parts = file.split('_')
                reaction_name = parts[3] if len(parts) >= 4 else file.split('.')[0]
                reaction_name = reaction_name.split('.')[0]
                data.append({
                    'reaction': reaction_name,
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

def process_input_csv(input_csv, reaction_col, ec_col):
    try:
        df = pd.read_csv(input_csv)
        return df[[reaction_col, ec_col]].copy()
    except Exception as e:
        print(f"Error reading input CSV: {e}")
        exit(1)

def main():
    parser = argparse.ArgumentParser(description="Process BridgIT results with any reference dataset")
    parser.add_argument("-b", "--bridgit_dir", required=True, help="Directory containing BridgIT zip results")
    parser.add_argument("-c", "--input_csv", required=True, help="Path to input CSV dataset")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file path")
    parser.add_argument("--bridgit_ec_col", default="reactionsA/ECA", help="Column in BridgIT txt files with EC predictions")
    parser.add_argument("--csv_reaction_col", required=True, help="Column in input CSV with reaction IDs")
    parser.add_argument("--csv_ec_col", required=True, help="Column in input CSV with EC numbers")
    args = parser.parse_args()

    bridgit_df = process_bridgit_results(args.bridgit_dir, args.bridgit_ec_col)
    input_df = process_input_csv(args.input_csv, args.csv_reaction_col, args.csv_ec_col)

    merged = pd.merge(
        bridgit_df, input_df,
        left_on="reaction", right_on=args.csv_reaction_col,
        how="inner"
    )

    final = merged[[args.csv_reaction_col, 'EC_number_predicted']].drop_duplicates()
    final.to_csv(args.output, index=False)
    print(f"Saved results to {args.output}")

if __name__ == "__main__":
    main()
