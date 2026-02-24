# -*- coding: utf-8 -*-
"""
Canonicalize reaction SMILES from a plain text file or CSV/TSV using RDKit.
"""

import os
import argparse
import pandas as pd
from rdkit import Chem


def canonicalize_mol_smiles(smiles: str) -> str | None:
    """Canonicalize a single molecule SMILES. Returns None if invalid."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:
        return None
    return None


def canonicalize_reaction_smiles(reaction_smiles: str) -> str | None:
    """
    Canonicalize a reaction SMILES (no atom mapping assumed).
    Returns standardized reaction SMILES: canonical(reactants)>>canonical(products)
    """
    try:
        reactants, products = reaction_smiles.strip().split(">>")
    except ValueError:
        print(f"Skipping invalid reaction format: {reaction_smiles}")
        return None

    reactant_cans = sorted(filter(None, [canonicalize_mol_smiles(s) for s in reactants.split('.') if s]))
    product_cans  = sorted(filter(None, [canonicalize_mol_smiles(s) for s in products.split('.')  if s]))

    if not reactant_cans and not product_cans:
        return None

    return '.'.join(reactant_cans) + ">>" + '.'.join(product_cans)


def detect_separator(path: str) -> str:
    """Infer separator from file extension."""
    return '\t' if path.endswith('.tsv') else ','


def process_tabular(
    input_path: str,
    output_path: str,
    smiles_col: str,
    output_col: str = "canonical_smiles",
    sep: str | None = None,
    keep_failed: bool = False,
) -> tuple[int, int]:
    """
    Canonicalize a SMILES column in a CSV/TSV file.
    Writes all original columns plus a new canonical SMILES column.
    Returns (n_success, n_failed).
    """
    sep = sep or detect_separator(input_path)
    df  = pd.read_csv(input_path, sep=sep)

    if smiles_col not in df.columns:
        available = ', '.join(df.columns.tolist())
        raise ValueError(f"Column '{smiles_col}' not found. Available columns: {available}")

    df[output_col] = df[smiles_col].astype(str).apply(canonicalize_reaction_smiles)

    failed  = df[output_col].isna().sum()
    success = len(df) - failed

    if not keep_failed:
        df = df.dropna(subset=[output_col])

    out_sep = sep if sep != ',' else ','
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    df.to_csv(output_path, sep=out_sep, index=False)

    return success, failed


def process_text_file(input_path: str, output_path: str) -> tuple[int, int]:
    """
    Canonicalize a plain text file where each line is a reaction SMILES.
    Returns (n_success, n_failed).
    """
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    success, failed = 0, 0

    with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            can_rxn = canonicalize_reaction_smiles(line)
            if can_rxn:
                f_out.write(can_rxn + '\n')
                success += 1
            else:
                print(f"Warning: could not canonicalize: {line}")
                failed += 1

    return success, failed


def process_input(
    input_path: str,
    output_path: str,
    smiles_col: str | None = None,
    output_col: str = "canonical_smiles",
    sep: str | None = None,
    keep_failed: bool = False,
    extension: str = ".txt",
) -> None:
    """
    Dispatch to the appropriate processor based on input type (file or directory).
    """
    if os.path.isdir(input_path):
        for root, _, files in os.walk(input_path):
            for fname in files:
                if not fname.endswith(extension):
                    continue
                in_file  = os.path.join(root, fname)
                out_file = os.path.join(output_path, os.path.relpath(in_file, input_path))
                print(f"Processing: {in_file}")
                success, failed = _process_single(in_file, out_file, smiles_col, output_col, sep, keep_failed)
                print(f"  Done: {success} canonicalized, {failed} failed -> {out_file}")

    elif os.path.isfile(input_path):
        success, failed = _process_single(input_path, output_path, smiles_col, output_col, sep, keep_failed)
        print(f"Done: {success} canonicalized, {failed} failed -> {output_path}")

    else:
        raise FileNotFoundError(f"Input path not found: {input_path}")


def _process_single(
    input_path: str,
    output_path: str,
    smiles_col: str | None,
    output_col: str,
    sep: str | None,
    keep_failed: bool,
) -> tuple[int, int]:
    """Route a single file to tabular or plain-text processor."""
    is_tabular = input_path.endswith(('.csv', '.tsv'))

    if is_tabular:
        if not smiles_col:
            raise ValueError(f"--smiles_col is required for tabular files (got: {input_path})")
        return process_tabular(input_path, output_path, smiles_col, output_col, sep, keep_failed)
    else:
        return process_text_file(input_path, output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Canonicalize reaction SMILES from a text file, CSV, or TSV using RDKit."
    )
    parser.add_argument("--input",        required=True,  help="Input file or directory.")
    parser.add_argument("--output",       required=True,  help="Output file or directory.")
    parser.add_argument("--smiles_col",   default=None,   help="Column name containing SMILES (required for CSV/TSV).")
    parser.add_argument("--output_col",   default="canonical_smiles", help="Output column name for canonical SMILES (default: canonical_smiles).")
    parser.add_argument("--sep",          default=None,   help="Delimiter for tabular files. Auto-detected from extension if not set.")
    parser.add_argument("--keep_failed",  action="store_true", help="Keep rows where canonicalization failed (written as empty).")
    parser.add_argument("--extension",    default=".txt", help="File extension to process when input is a directory (default: .txt).")
    args = parser.parse_args()

    process_input(
        input_path   = args.input,
        output_path  = args.output,
        smiles_col   = args.smiles_col,
        output_col   = args.output_col,
        sep          = args.sep,
        keep_failed  = args.keep_failed,
        extension    = args.extension,
    )


if __name__ == "__main__":
    main()