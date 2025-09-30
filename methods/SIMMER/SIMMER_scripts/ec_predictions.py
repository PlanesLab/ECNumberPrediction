import os
import pandas as pd
import argparse

def process_file(file_path):
    """
    Process a single TSV file and return all valid EC numbers that have a sub-subclass prediction.
    A valid EC number is one that, when split by '.', has at least three components and the third component is non-empty and not a placeholder.
    Returns a semicolon-separated string of valid EC numbers, or None if none are found.
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None

    if 'EC' not in df.columns:
        print(f"Column 'EC' not found in file {file_path}.")
        return None

    # Ensure the EC column is of type string
    df['EC'] = df['EC'].astype(str)
    
    valid_ecs = []
    for ec in df['EC']:
        try:
            ec_split = ec.split('.')
            if len(ec_split) >= 3 and ec_split[2] not in ['', '-', 'None']:
                valid_ecs.append(ec)
        except Exception as e:
            print(f"Error processing EC number '{ec}' in file {file_path}: {e}")
    return ';'.join(valid_ecs) if valid_ecs else None


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process SIMMER EC predictions.")
    parser.add_argument('--input_dir', required=True, help="Directory containing input TSV files")
    parser.add_argument('--output_file', required=True, help="Path to save the output CSV file")
    args = parser.parse_args()

    # Get all TSV files in the input directory
    try:
        files = [f for f in os.listdir(args.input_dir) if f.endswith('.tsv')]
    except Exception as e:
        print(f"Error listing files in {args.input_dir}: {e}")
        return

    data = []
    for file in files:
        reaction_name = file.split('_')[0].strip()
        file_path = os.path.join(args.input_dir, file)
        simmer_pred = process_file(file_path)
        data.append({'Reaction name': reaction_name, 'SIMMER pred': simmer_pred})
    
    df_predictions = pd.DataFrame(data)
    
    try:
        df_predictions.to_csv(args.output_file, index=False)
        print(f"Data saved to {args.output_file}")
    except Exception as e:
        print(f"Error saving file {args.output_file}: {e}")


if __name__ == "__main__":
    main()
