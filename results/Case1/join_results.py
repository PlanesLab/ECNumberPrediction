import pandas as pd
import glob
import os
from functools import reduce

# -------------------------------
# Part 1: Merge Query_Results CSVs
# -------------------------------

# Define the folder path and file pattern for CSV files in Query_Results
folder_path = "results/Case1/results"
file_pattern = os.path.join(folder_path, "*.csv")

# Get a list of all CSV files in the folder
files = glob.glob(file_pattern)

if not files:
    print(f"No CSV files found in {folder_path}")
    exit()

# List to hold each DataFrame from the Query_Results folder
dataframes = []

for file in files:
    # Read the CSV file (adjust 'sep' if needed)
    df = pd.read_csv(file, sep=',')
    
    # Rename the first column to 'reaction_id' regardless of its original name
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "reaction_id"})
    
    # Create a prefix based on the filename (without extension)
    prefix = os.path.splitext(os.path.basename(file))[0] + "_"
    
    # Rename all columns except the join key ('reaction_id') by prefixing them with the filename
    new_columns = {col: prefix + col for col in df.columns if col != "reaction_id"}
    df = df.rename(columns=new_columns)
    
    dataframes.append(df)


merged_df = reduce(lambda left, right: pd.merge(left, right, on="reaction_id", how="outer"), dataframes)

# -------------------------------
# Part 2: Clean the Merged DataFrame
# -------------------------------

def clean_value(x):
    # If the value is already missing (NaN), return None
    if pd.isna(x):
        return None
    # If the value is a string, check for unwanted patterns
    if isinstance(x, str):
        stripped = x.strip()
        if stripped in ["", "No Significant EC", "No EC Prediction" ,"No|EC|Prediction", "nan", "nan|", "|nan"]:
            return None
    return x

# Apply the cleaning function to every cell in the merged DataFrame
merged_df = merged_df.applymap(clean_value)

# -------------------------------
# Part 3: Merge with the KEGG CSV
# -------------------------------

# Define the KEGG CSV file path
kegg_file = "data/KEGG/kegg_reactions_current_test.csv"

# Read the KEGG CSV file (adjust 'sep' if needed)
kegg_df = pd.read_csv(kegg_file)

# Rename its first column to "reaction_id"
kegg_df = kegg_df.rename(columns={kegg_df.columns[0]: "reaction_id"})

print(kegg_df.columns)

# Keep only the "reaction_id" and "EC Number" columns
if "EC Number" in kegg_df.columns:
    kegg_df = kegg_df[["reaction_id", "EC Number"]]
else:
    print("Column 'EC Number' not found in the KEGG file.")
    exit()

# Rename the "EC Number" column to "KEGG"
kegg_df = kegg_df.rename(columns={"Reaction_EC_number": "KEGG"})

# (Optional) Clean the KEGG dataframe as well using the same cleaning function
kegg_df = kegg_df.applymap(clean_value)

# Merge the KEGG dataframe with the merged_df from Query_Results on "reaction_id"
merged_df = pd.merge(merged_df, kegg_df, on="reaction_id", how="outer")

# -------------------------------
# Part 4: Save the Final Merged DataFrame
# -------------------------------

# Define the output file path
output_file = os.path.join(folder_path, "results/Case1/merged_output.csv")

# Save the final merged and cleaned DataFrame to a CSV file
merged_df.to_csv(output_file, index=False)

print(f"Merged and cleaned file saved as '{output_file}'.")