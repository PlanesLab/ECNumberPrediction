"""
Script: query_selenzyme.py
Author: Josefina Arcagni
Date: 2025-09-09
Description: Query the Selenzyme server with reactions from a CSV file and save results.
"""

import argparse
import csv
import os

import pandas as pd
import requests
from bs4 import BeautifulSoup


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Query Selenzyme with reactions from a CSV file.")
    parser.add_argument("--server_url", type=str, default="http://localhost:32784",
                        help="Base URL for Selenzyme server (default: http://localhost:32784).")
    parser.add_argument("--csv_file", type=str, required=True,
                        help="Path to input CSV file.")
    parser.add_argument("--results_folder", type=str, required=True,
                        help="Folder to save per-reaction query results.")
    parser.add_argument("--reaction_name_column", type=str, default="drug",
                        help="Column name for the reaction identifier (default: drug).")
    parser.add_argument("--reaction_smiles_column", type=str, default="reaction_smiles",
                        help="Column name for the reaction SMILES (default: reaction_smiles).")
    return parser.parse_args()


def query_reaction(
    session: requests.Session,
    base_url: str,
    reaction_id: str,
    smiles: str,
    output_file: str,
) -> None:
    smarts_data = {"smarts": smiles, "rdb": "ec", "rxnid": ""}
    display_response = session.post(f"{base_url}/display", data=smarts_data)
    if display_response.status_code != 200:
        print(f"Step 1 failed for {reaction_id} (status {display_response.status_code})")
        return

    results_data = {"targets": "200", "noMSA": "on", "host": "83333", "finger": "Morgan"}
    results_response = session.post(f"{base_url}/results", data=results_data)
    if results_response.status_code not in [200, 302]:
        print(f"Step 2 failed for {reaction_id} (status {results_response.status_code})")
        return

    soup = BeautifulSoup(results_response.text, "html.parser")
    table = soup.find("table")
    if not table:
        print(f"No results table for reaction {reaction_id}.")
        return

    rows_html = table.find_all("tr")
    table_data = [
        [cell.get_text(strip=True) for cell in row.find_all(["th", "td"])]
        for row in rows_html
    ]
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(table_data)
    print(f"Saved results for {reaction_id} → {output_file}")


def main() -> None:
    args = parse_args()
    os.makedirs(args.results_folder, exist_ok=True)

    df = pd.read_csv(args.csv_file)

    for _, row in df.iterrows():
        reaction_id = str(row[args.reaction_name_column])
        output_file = os.path.join(args.results_folder, f"{reaction_id}.csv")

        if os.path.exists(output_file):
            print(f"Skipping {reaction_id} – already processed.")
            continue

        smiles = str(row[args.reaction_smiles_column])
        print(f"\nProcessing reaction {reaction_id}...")

        session = requests.Session()
        session.get(f"{args.server_url}/")
        query_reaction(session, args.server_url, reaction_id, smiles, output_file)

    print("\nDone.")


if __name__ == "__main__":
    main()
