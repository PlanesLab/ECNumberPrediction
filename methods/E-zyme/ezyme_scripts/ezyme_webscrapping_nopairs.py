import os
import time
import itertools
import argparse
import tempfile
import requests
import pandas as pd
from bs4 import BeautifulSoup
from rdkit import Chem


EZYME_URL = "https://www.genome.jp/tools-bin/e-zyme-ko"


def split_reaction_smiles(reaction_smiles: str):
    reactants_part, products_part = reaction_smiles.split(">>")
    reactants = [r.strip() for r in reactants_part.split(".") if r.strip()]
    products = [p.strip() for p in products_part.split(".") if p.strip()]
    return reactants, products


def smiles_to_molfile(smiles: str, path: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    Chem.MolToMolFile(mol, path)


def run_ezyme_with_molfiles(reactant_smiles, product_smiles):
    with tempfile.TemporaryDirectory() as tmpdir:
        mol1_path = os.path.join(tmpdir, "reactant.mol")
        mol2_path = os.path.join(tmpdir, "product.mol")

        smiles_to_molfile(reactant_smiles, mol1_path)
        smiles_to_molfile(product_smiles, mol2_path)

        # STEP 1: view
        with open(mol1_path, "rb") as f1, open(mol2_path, "rb") as f2:
            r1 = requests.post(
                EZYME_URL,
                files={
                    "file1": ("reactant.mol", f1, "chemical/x-mdl-molfile"),
                    "file2": ("product.mol", f2, "chemical/x-mdl-molfile"),
                },
                data={"mode": "view"},
                timeout=60,
            )

        r1.raise_for_status()
        soup1 = BeautifulSoup(r1.text, "html.parser")

        id_input = soup1.find("input", {"name": "id"})
        if not id_input:
            raise RuntimeError("Failed to obtain E-zyme session ID")

        ez_id = id_input["value"]

        # STEP 2: compute
        with open(mol1_path, "rb") as f1, open(mol2_path, "rb") as f2:
            r2 = requests.post(
                EZYME_URL,
                files={
                    "file1": ("reactant.mol", f1, "chemical/x-mdl-molfile"),
                    "file2": ("product.mol", f2, "chemical/x-mdl-molfile"),
                },
                data={"mode": "compute", "id": ez_id},
                timeout=60,
            )

        r2.raise_for_status()
        return BeautifulSoup(r2.text, "html.parser")


def extract_tables(soup):
    results = {}

    # E-zyme2
    div2 = soup.find("div", id="ref_rp_img")
    if div2 and div2.find("table"):
        rows = div2.find("table").find_all("tr")[1:]
        data = [[c.get_text(" ", strip=True) for c in r.find_all("td")] for r in rows]
        results["EZYME2"] = pd.DataFrame(
            data, columns=["RPAIR", "Score", "EC", "KO"]
        )

    # E-zyme1
    div1 = soup.find("div", id="ez1")
    if div1 and div1.find("table"):
        rows = div1.find("table").find_all("tr")[1:]
        data = [[c.get_text(" ", strip=True) for c in r.find_all("td")] for r in rows]
        results["EZYME1"] = pd.DataFrame(
            data, columns=["EC", "Weighted Score", "Observed Freq", "Reactions"]
        )

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_tsv", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("--sleep", type=float, default=2.0)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = pd.read_csv(args.input_tsv, sep="\t")

    for _, row in df.iterrows():
        reaction_id = str(row["ID"])
        reaction_smiles = row["REACTION_SMILES"]

        try:
            reactants, products = split_reaction_smiles(reaction_smiles)
        except Exception as e:
            print(f"❌ {reaction_id}: {e}")
            continue

        pairs = list(itertools.product(reactants, products))
        print(f"🔬 Reaction {reaction_id}: {len(pairs)} pairs")

        reaction_dir = os.path.join(args.output_dir, reaction_id)
        os.makedirs(reaction_dir, exist_ok=True)

        for i, (r_smiles, p_smiles) in enumerate(pairs, start=1):
            pair_dir = os.path.join(reaction_dir, f"pair_{i}")
            os.makedirs(pair_dir, exist_ok=True)

            try:
                soup = run_ezyme_with_molfiles(r_smiles, p_smiles)
                tables = extract_tables(soup)

                for name, df_out in tables.items():
                    df_out.to_csv(os.path.join(pair_dir, f"{name}.csv"), index=False)

                print(f"✅ {reaction_id} pair_{i}")

            except Exception as e:
                print(f"❌ {reaction_id} pair_{i}: {e}")

            time.sleep(args.sleep)

    print("✅ All done")


if __name__ == "__main__":
    main()
