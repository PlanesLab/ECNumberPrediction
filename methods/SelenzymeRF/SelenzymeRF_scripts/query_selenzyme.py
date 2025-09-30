"""
Script: query_selenzyme.py
Author: Josefina Arcagni
Date: 2025-09-09
Description: This script queries the Selenzyme server with reactions from a CSV file and saves the results.
"""

import os
import csv
import pandas as pd
import requests
from bs4 import BeautifulSoup
import argparse

parser = argparse.ArgumentParser(description="Query Selenzyme with reactions from a CSV file.")
parser.add_argument("--base_url", type=str, default="http://localhost:32784", help="Base URL for Selenzyme server")
parser.add_argument("--csv_file", type=str, help="Path to input CSV file")
parser.add_argument("--results_folder", type=str, help="Folder to save query results")
parser.add_argument("--reaction_name_column", type=str, default="drug", help="Column name for reaction name")
parser.add_argument("--reaction_smiles_column", type=str, default="reaction_smiles", help="Column name for reaction SMILES")
args = parser.parse_args()

BASE_URL = args.base_url
CSV_FILE = args.csv_file
RESULTS_FOLDER = args.results_folder
REACTION_NAME_COL = args.reaction_name_column
REACTION_SMILES_COL = args.reaction_smiles_column

os.makedirs(RESULTS_FOLDER, exist_ok=True)

# Read the CSV file using pandas (assumes headers are present)
df = pd.read_csv(CSV_FILE)

# Loop over each reaction (each row) in the CSV
for index, row in df.iterrows():
    reaction_id = row[REACTION_NAME_COL]
    # Construct output path and skip if it already exists
    output_file = os.path.join(RESULTS_FOLDER, f"{reaction_id}.csv")
    if os.path.exists(output_file):
        print(f"Skipping reaction {reaction_id}: results already exist.")
        continue

    isomeric_smiles = row[REACTION_SMILES_COL]
    print(f"\nProcessing reaction {reaction_id}...")

    session = requests.Session()

    home_response = session.get(f"{BASE_URL}/")
    soup_home = BeautifulSoup(home_response.text, "html.parser")
    # CSRF token extraction (if needed) is commented out

    smarts_data = {
        "smarts": isomeric_smiles,
        "rdb": "ec",
        "rxnid": ""
    }
    display_url = f"{BASE_URL}/display"
    display_response = session.post(display_url, data=smarts_data)
    if display_response.status_code == 200:
        print("Step 1: Reaction input submitted successfully!")
    else:
        print(f"Step 1: Failed to submit reaction input. (Status Code: {display_response.status_code})")
        continue

    results_data = {
        "targets": "200",
        "noMSA": "on",
        "host": "83333",
        "finger": "Morgan"
    }
    results_url = f"{BASE_URL}/results"
    results_response = session.post(results_url, data=results_data)
    if results_response.status_code in [200, 302]:
        print("Step 2: Results submitted successfully!")
    else:
        print(f"Step 2: Failed to submit results. (Status Code: {results_response.status_code})")
        continue

    soup_results = BeautifulSoup(results_response.text, "html.parser")
    table = soup_results.find("table")
    if not table:
        print(f"Step 3: No table found in the response for reaction {reaction_id}.")
        continue

    rows_html = table.find_all("tr")
    table_data = []
    for row_html in rows_html:
        cells = row_html.find_all(["th", "td"])
        row_data = [cell.get_text(strip=True) for cell in cells]
        table_data.append(row_data)

    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerows(table_data)
    print(f"Step 4: Table data saved to {output_file}")
