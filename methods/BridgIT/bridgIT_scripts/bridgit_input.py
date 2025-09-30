# ==============================================================================
# Author: Josefina Arcagni
# Email: jarcagniriv@unav.es
# Date: 2025-02-07
# Script: bridgit_input.py
# Description: Prepares KEGG data for input in the BridgIT LCSB server.
# ==============================================================================
import os
import csv
import requests
import re
import time
import argparse

parser = argparse.ArgumentParser(description="Prepare data for BridgIT input.")
parser.add_argument('--input_file', required=True, help='Path to input CSV file')
parser.add_argument('--output_file', required=True, help='Path to output formatted reactions file')
parser.add_argument('--molfile_folder', required=True, help='Folder to save retrieved molfiles')
parser.add_argument('--equation_column', default='Equation',
                    help='Name of the column containing reaction equations (default: Equation)')
args = parser.parse_args()

input_file = args.input_file
output_file = args.output_file
molfile_folder = args.molfile_folder
equation_column = args.equation_column

os.makedirs(molfile_folder, exist_ok=True)
retrieved_compounds = set()
delay_time = 1  

# ----------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------
def format_equation(equation: str) -> str:
    formatted_equation = []
    for part in equation.split():
        main_part = part.split('(')[0]
        if main_part.isdigit():
            formatted_equation.append(f"({main_part})")
        elif re.match(r"^[CG]\d{5}$", main_part):
            formatted_equation.append(part)
        elif part in {"+", "<=>", "=>"}:
            formatted_equation.append(part)
        else:
            print(f"Unexpected element found: {part}")
    return ' '.join(formatted_equation)

def extract_kegg_compounds(equation: str) -> set:
    return set(re.findall(r'\b[CG]\d{5}\b', equation))

def retrieve_molfile(compound_id: str) -> None:
    if compound_id in retrieved_compounds:
        return
    url = f"http://rest.kegg.jp/get/{compound_id}/mol"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            molfile_path = os.path.join(molfile_folder, f"{compound_id}.mol")
            with open(molfile_path, mode='w') as molfile:
                molfile.write(response.text)
            retrieved_compounds.add(compound_id)
            print(f"Saved molfile for {compound_id}")
        else:
            print(f"Failed to retrieve molfile for {compound_id} (HTTP {response.status_code})")
    except requests.RequestException as e:
        print(f"Error retrieving molfile for {compound_id}: {e}")
    time.sleep(delay_time)

# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    with open(input_file, mode='r') as infile:
        reader = csv.DictReader(infile)
        entries = []
        for index, row in enumerate(reader, start=1):
            if equation_column not in row:
                raise KeyError(f"Column '{equation_column}' not found in input file")
            equation = row[equation_column]
            formatted_equation = format_equation(equation)
            compounds = extract_kegg_compounds(equation)
            for compound in compounds:
                retrieve_molfile(compound)
            entries.append(f"{index};;{formatted_equation};")

    with open(output_file, mode='w') as outfile:
        outfile.write("ENTRY;KEGG;EQUATION;OPERATORS\n")
        outfile.write("\n".join(entries))

    print(f"[DONE] File saved to {output_file}")

if __name__ == "__main__":
    main()
