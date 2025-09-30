import subprocess
import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser(description="Query Theia with reaction queries and IDs.")
parser.add_argument("--query_file", required=True, help="Path to the file containing queries.")
parser.add_argument("--reaction_ids_file", required=True, help="Path to the CSV file containing reaction IDs.")
parser.add_argument("--output_file", required=True, help="Path to the output CSV file.")
parser.add_argument("--reaction_id_column", required=True, help="Name of the column containing reaction IDs.")
parser.add_argument("--model", default="ecreact.ec123", help="Theia model to use (default: ecreact.ec123).")

args = parser.parse_args()

query_file = args.query_file
reaction_ids_file = args.reaction_ids_file
output_file = args.output_file
reaction_id_column = args.reaction_id_column
model = args.model

with open(query_file, 'r') as f:
    queries = [line.strip() for line in f if line.strip()]

reaction_df = pd.read_csv(reaction_ids_file)
if reaction_id_column not in reaction_df.columns:
    raise ValueError(f"The reaction ids file must contain a '{reaction_id_column}' column.")

reaction_ids = reaction_df[reaction_id_column].tolist()

if len(queries) != len(reaction_ids):
    raise ValueError("The number of queries and reaction IDs do not match.")

predictions = []
for query in queries:
    command = ["theia-cli", model, query, "--probs"]
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        prediction = result.stdout.strip()
    except subprocess.CalledProcessError as e:
        prediction = f"Error: {e.stderr.strip()}"
    predictions.append(prediction)

result_df = pd.DataFrame({
    reaction_id_column: reaction_ids,
    "Prediction": predictions
})

result_df.to_csv(output_file, index=False)
print(f"Results saved to {output_file}")
