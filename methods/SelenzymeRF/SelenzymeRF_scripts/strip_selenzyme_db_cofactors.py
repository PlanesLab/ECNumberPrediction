"""
Build a cofactor-stripped copy of the SelenzymeRF reference database
(data_2023/), for the oxidoreductase no-cofactor ablation.

Strips reac_smi.csv via SIMMER's strip_cofactor_components. FP_Morg.npz
(per-compound) needs no changes. FP_MorgRF.npz rows for (cofactor compound,
modified reaction) pairs are dropped rather than recomputed -- an
approximation, since reacting-fragment fingerprints were never re-derived
via RXNMapper.

Cofactor identification uses chem_prop.tsv (MNXM_ID -> SMILES) from MNXref
4.5 (2025-08), since the bundled 2023 snapshot lacks it.

Usage:
    python strip_selenzyme_db_cofactors.py \
        --base_db_dir methods/SelenzymeRF/SelenzymeRF_code/data_2023 \
        --output_dir methods/SelenzymeRF/SelenzymeRF_code/data_2023_nocofactor
"""

import argparse
import os
import shutil
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "SIMMER", "SIMMER_scripts"))
from cofactors import COFACTOR_KEYS, _skeleton_key, strip_cofactor_components  # noqa: E402


def get_referenced_compound_ids(reac_prop_path: str) -> set:
    ids = set()
    with open(reac_prop_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 2:
                continue
            left, right = row[1].split(" = ")
            for side in (left, right):
                for term in side.split(" + "):
                    cid = term.strip().split(" ")[-1].split("@")[0]
                    ids.add(cid)
    return ids


def lookup_smiles(chem_prop_path: str, wanted_ids: set) -> dict:
    id_to_smiles = {}
    with open(chem_prop_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            tab_idx = line.find("\t")
            if tab_idx == -1:
                continue
            cid = line[:tab_idx]
            if cid not in wanted_ids:
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 9:
                continue
            smiles = row[8].strip()
            if smiles:
                id_to_smiles[cid] = smiles
    return id_to_smiles


def identify_cofactor_ids(id_to_smiles: dict) -> set:
    cofactor_ids = set()
    for cid, smi in id_to_smiles.items():
        key = _skeleton_key(smi)
        if key is not None and key in COFACTOR_KEYS:
            cofactor_ids.add(cid)
    return cofactor_ids


def strip_reac_prop(reac_prop_path: str, out_path: str, cofactor_ids: set):
    n_total = n_modified = n_dropped = 0
    modified_reactions = set()
    with open(reac_prop_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            row = line.rstrip("\n").split("\t")
            n_total += 1
            rid = row[0]
            left, right = row[1].split(" = ")

            def strip_side(side):
                kept_terms = []
                for term in side.split(" + "):
                    parts = term.strip().split(" ")
                    cid = parts[-1].split("@")[0]
                    if cid not in cofactor_ids:
                        kept_terms.append(term.strip())
                return kept_terms

            left_kept = strip_side(left)
            right_kept = strip_side(right)
            was_modified = (len(left_kept) != len(left.split(" + "))) or \
                            (len(right_kept) != len(right.split(" + ")))

            if not left_kept or not right_kept:
                n_dropped += 1
                continue
            if was_modified:
                n_modified += 1
                modified_reactions.add(rid)
                row[1] = " + ".join(left_kept) + " = " + " + ".join(right_kept)
            fout.write("\t".join(row) + "\n")
    print(f"reac_prop.tsv: total={n_total} modified={n_modified} dropped={n_dropped}")
    return modified_reactions


def strip_reac_smi(reac_smi_path: str, out_path: str):
    df = pd.read_csv(reac_smi_path, dtype=str, keep_default_na=False)
    kept_rows = []
    n_modified = n_dropped = 0
    for _, row in df.iterrows():
        try:
            left, right = row["SMILES"].split(">>")
        except ValueError:
            kept_rows.append(row)
            continue
        left_kept, left_dropped = strip_cofactor_components(left.strip())
        right_kept, right_dropped = strip_cofactor_components(right.strip())
        if not left_kept or not right_kept:
            n_dropped += 1
            continue
        if left_dropped or right_dropped:
            n_modified += 1
        new_row = row.copy()
        new_row["SMILES"] = f"{left_kept}>>{right_kept}"
        kept_rows.append(new_row)
    pd.DataFrame(kept_rows).to_csv(out_path, index=False)
    print(f"reac_smi.csv: total={len(df)} modified={n_modified} dropped={n_dropped} kept={len(kept_rows)}")


def filter_fp_morgrf(fp_morgrf_path: str, out_path: str, cofactor_ids: set, modified_reactions: set):
    data = np.load(fp_morgrf_path, allow_pickle=True)
    x, y, z, d = data["x"], data["y"], data["z"], data["d"]
    keep_mask = ~(np.isin(y, list(cofactor_ids)) & np.isin(z, list(modified_reactions)))
    n_total = len(y)
    n_dropped = n_total - keep_mask.sum()
    np.savez_compressed(out_path, x=x[keep_mask], y=y[keep_mask], z=z[keep_mask], d=d[keep_mask])
    print(f"FP_MorgRF.npz: total_rows={n_total} dropped_cofactor_rows={n_dropped}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_db_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--chem_prop", default=None,
                         help="Path to chem_prop.tsv (default: <base_db_dir>/chem_prop.tsv)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    chem_prop_path = args.chem_prop or os.path.join(args.base_db_dir, "chem_prop.tsv")

    print("1. Identifying referenced compound ids...")
    wanted_ids = get_referenced_compound_ids(os.path.join(args.base_db_dir, "reac_prop.tsv"))
    print(f"   {len(wanted_ids)} unique compound ids referenced in reac_prop.tsv")

    print("2. Looking up SMILES from chem_prop.tsv...")
    id_to_smiles = lookup_smiles(chem_prop_path, wanted_ids)
    print(f"   resolved {len(id_to_smiles)} / {len(wanted_ids)}")

    print("3. Identifying cofactor compound ids...")
    cofactor_ids = identify_cofactor_ids(id_to_smiles)
    print(f"   {len(cofactor_ids)} cofactor compound ids identified")
    with open(os.path.join(args.output_dir, "cofactor_ids.txt"), "w") as f:
        for cid in sorted(cofactor_ids):
            f.write(f"{cid}\t{id_to_smiles[cid]}\n")

    print("4. Stripping reac_prop.tsv...")
    modified_reactions = strip_reac_prop(
        os.path.join(args.base_db_dir, "reac_prop.tsv"),
        os.path.join(args.output_dir, "reac_prop.tsv"),
        cofactor_ids,
    )

    print("5. Stripping reac_smi.csv...")
    strip_reac_smi(
        os.path.join(args.base_db_dir, "reac_smi.csv"),
        os.path.join(args.output_dir, "reac_smi.csv"),
    )

    print("6. Filtering FP_MorgRF.npz...")
    filter_fp_morgrf(
        os.path.join(args.base_db_dir, "FP_MorgRF.npz"),
        os.path.join(args.output_dir, "FP_MorgRF.npz"),
        cofactor_ids, modified_reactions,
    )

    print("7. Passing through unchanged files (symlink)...")
    for filename in ["FP_Morg.npz", "reac_seqs.tsv", "reac_xref.tsv", "seq_org.tsv",
                      "org_lineage.csv", "rxn_consensus_20160612.txt"]:
        src = os.path.join(args.base_db_dir, filename)
        dst = os.path.join(args.output_dir, filename)
        if os.path.exists(src):
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(os.path.abspath(src), dst)
            print(f"   {filename}: symlinked")
        else:
            print(f"   WARNING: {src} not found, skipping")

    print(f"\nCofactor-stripped base DB written to {args.output_dir}")


if __name__ == "__main__":
    main()
