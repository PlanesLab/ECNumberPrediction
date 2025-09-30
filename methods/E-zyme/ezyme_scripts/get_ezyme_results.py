import os
import pandas as pd
import argparse

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description="Collect E-zyme predictions from multiple folders.")
    parser.add_argument("--input_dir", required=True, help="Path to the root directory containing reaction folders")
    parser.add_argument("--output_file", required=True, help="Path to save the final CSV file")
    args = parser.parse_args()

    root_dir = args.input_dir
    csv_data = []

    for reaction_folder in os.listdir(root_dir):
        reaction_folder_path = os.path.join(root_dir, reaction_folder)

        if os.path.isdir(reaction_folder_path):
            reaction_files = os.listdir(reaction_folder_path)

            ezyme1_predictions = []
            ezyme2_predictions = []

            for file in reaction_files:
                if '_E-zyme1.csv' in file:
                    ezyme1_file = os.path.join(reaction_folder_path, file)
                    try:
                        ezyme1_df = pd.read_csv(ezyme1_file)
                        preds = ezyme1_df.iloc[:, 0].astype(str).str.replace(" ", "|").tolist()
                        ezyme1_predictions.extend(preds)
                    except Exception as e:
                        print(f"Error reading E-zyme1 file in {reaction_folder}: {e}")
                elif '_E-zyme2.csv' in file:
                    ezyme2_file = os.path.join(reaction_folder_path, file)
                    try:
                        ezyme2_df = pd.read_csv(ezyme2_file)
                        preds = ezyme2_df.iloc[:, 2].astype(str).str.replace(" ", "|").tolist()
                        ezyme2_predictions.extend(preds)
                    except Exception as e:
                        print(f"Error reading E-zyme2 file in {reaction_folder}: {e}")

            ezyme1_result = ';'.join(ezyme1_predictions) if ezyme1_predictions else "No EC Prediction"
            ezyme2_result = ';'.join(ezyme2_predictions) if ezyme2_predictions else "No EC Prediction"

            csv_data.append([reaction_folder, ezyme1_result, ezyme2_result])

    # Save final CSV
    final_df = pd.DataFrame(csv_data, columns=['Reaction Name', 'E-zyme1', 'E-zyme2'])
    final_df.to_csv(args.output_file, index=False)

    print(f"CSV file created successfully: {args.output_file}")


if __name__ == "__main__":
    main()
