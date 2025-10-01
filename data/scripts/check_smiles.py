# -*- coding: utf-8 -*-
"""
Canonicalize all reaction SMILES in a file
"""

from rdkit import Chem
import os

def canonicalize_mol_smiles(smiles):
    """
    Canonicalize a single molecule SMILES.
    Returns None if the SMILES is invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Remove kekulization to preserve aromatics/stereochem if needed
            return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:
        return None
    return None

def canonicalize_reaction_smiles(reaction_smiles):
    """
    Canonicalize a reaction SMILES (no atom mapping assumed).
    Returns standardized reaction SMILES in the form:
        canonical(reactant1.reactant2)>>canonical(product1.product2)
    """
    try:
        reactants, products = reaction_smiles.strip().split(">>")
    except ValueError:
        print(f"Skipping invalid reaction format: {reaction_smiles}")
        return None

    reactant_mols = reactants.split('.') if reactants else []
    product_mols = products.split('.') if products else []

    reactant_cans = [canonicalize_mol_smiles(smi) for smi in reactant_mols]
    product_cans = [canonicalize_mol_smiles(smi) for smi in product_mols]

    # Filter out invalid molecules
    reactant_cans = [s for s in reactant_cans if s]
    product_cans = [s for s in product_cans if s]

    reactant_cans.sort()
    product_cans.sort()

    return '.'.join(reactant_cans) + ">>" + '.'.join(product_cans)

def process_file(input_path, output_path):
    """
    Read a file of reaction SMILES, canonicalize them, and write to a new file.
    """
    with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            can_rxn = canonicalize_reaction_smiles(line)
            if can_rxn:
                f_out.write(can_rxn + '\n')
            else:
                print(f"Warning: could not canonicalize line:\n{line}")

# Paths
input_file = "/scratch/jarcagniriv/CaseStudy/drugs/reaction_smiles.txt"
output_file = "/scratch/jarcagniriv/CaseStudy/drugs/reaction_smiles_can.txt"

# Run
process_file(input_file, output_file)
print(f"Canonicalized reactions written to: {output_file}")
