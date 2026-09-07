"""
Split a bridgit_input.py-generated systemfile.txt + molfiles/ folder into small
per-batch ZIPs for upload to https://lcsb-databases.epfl.ch/Bridgit (which has a
practical per-submission size/reaction-count limit -- the original Case 1 KEGG-1.8K
input used ~37 reactions/batch across 50 ZIPs). Generalizes the KEGG/Rhea-hardcoded
one-off methods/BridgIT/bridgIT_scripts/bridgit_reduced_input.py into a reusable,
argparse'd script so it can be pointed at any seed's input dir.

Each output ZIP <output_dir>/reducedinput<i>.zip contains:
  reduced_systemfile.txt   (header + this batch's reaction rows)
  molfiles/*.mol           (only the compounds referenced by this batch)
"""

import argparse
import math
import os
import re
import shutil
import zipfile


def clean_equation(equation: str) -> str:
    equation = re.sub(r"\(n[+-]?\d*\)", "", equation)
    equation = re.sub(r"\(m[+-]?\d*\)", "", equation)
    return equation.replace(" ", "")


def main() -> None:
    parser = argparse.ArgumentParser(description="Split a BridgIT systemfile+molfiles into per-batch ZIPs.")
    parser.add_argument("--systemfile", required=True)
    parser.add_argument("--molfiles_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--batch_size", type=int, default=10, help="Reactions per ZIP (default 10, matching the ~37-reactions/50-batches ratio used for the original 1866-reaction KEGG-1.8K input).")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.systemfile) as f:
        lines = f.readlines()
    # bridgit_input.py writes a fixed 4-line preamble (COMPOUNDS/ENTRY/reactionsS/
    # ENTRY;KEGG;EQUATION;OPERATORS) before the actual reaction rows -- find the CSV header
    # line rather than assuming line count, so this still works if the preamble ever changes.
    header_end = next(i for i, line in enumerate(lines) if line.startswith("ENTRY;KEGG"))
    header = lines[0:header_end + 1]
    reactions = lines[header_end + 1:]
    num_reactions = len(reactions)
    num_batches = math.ceil(num_reactions / args.batch_size)

    for i in range(num_batches):
        start_idx = i * args.batch_size
        end_idx = min((i + 1) * args.batch_size, num_reactions)
        selected_reactions = []
        molecule_ids = set()

        for reaction in reactions[start_idx:end_idx]:
            parts = reaction.strip().split(";")
            if len(parts) > 2:
                parts[2] = clean_equation(parts[2])
            selected_reactions.append(";".join(parts) + "\n")

            for token in parts[2].split("+"):
                token = token.split("<=>")[0]
                token = token.replace("(", "").replace(")", "").strip()
                if token.isalnum():
                    molecule_ids.add(token)

        split_folder = os.path.join(args.output_dir, f"reducedinput{i + 1}")
        split_molfiles_folder = os.path.join(split_folder, "molfiles")
        os.makedirs(split_molfiles_folder, exist_ok=True)

        n_found = 0
        for mol_id in molecule_ids:
            source_path = os.path.join(args.molfiles_dir, f"{mol_id}.mol")
            if os.path.exists(source_path):
                shutil.copy(source_path, os.path.join(split_molfiles_folder, f"{mol_id}.mol"))
                n_found += 1

        output_file_path = os.path.join(split_folder, "reduced_systemfile.txt")
        with open(output_file_path, "w") as outfile:
            outfile.writelines(header)
            outfile.writelines(selected_reactions)

        zip_file_path = os.path.join(args.output_dir, f"reducedinput{i + 1}.zip")
        with zipfile.ZipFile(zip_file_path, "w", zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(split_folder):
                for file in files:
                    file_path = os.path.join(root, file)
                    zipf.write(file_path, os.path.relpath(file_path, split_folder))

        shutil.rmtree(split_folder)
        print(f"Batch {i + 1}/{num_batches}: {end_idx - start_idx} reactions, "
              f"{n_found}/{len(molecule_ids)} molfiles found -> {zip_file_path}")


if __name__ == "__main__":
    main()
