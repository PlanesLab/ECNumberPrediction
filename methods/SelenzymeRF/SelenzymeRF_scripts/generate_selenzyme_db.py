"""
Script: generate_selenzyme_db.py
Author: Josefina Arcagni
Date: 2025-09-09 (rewritten 2026-08-15)
Description:
    Build a filtered copy of the SelenzymeRF reference database (data_2023/)
    that excludes reactions in a given seed's test split, to prevent leakage
    when querying the server against that split.

    Filters reac_prop.tsv, reac_seqs.tsv, reac_smi.csv (keyed by MNXR reaction
    id; everything else the server uses derives from these). FP_Morg.npz and
    FP_MorgRF.npz are passed through unchanged: FP_Morg.npz is per-compound,
    not per-reaction, and FP_MorgRF.npz rows for excluded reactions are simply
    never looked up.

Usage:
    python generate_selenzyme_db.py \
        --base_db_dir /path/to/unzipped/data_2023 \
        --test_tsv data/Splits-DBs/MetaNetX/seed_splits/seed0/test.tsv \
        --id_scheme metanetx \
        --output_dir methods/SelenzymeRF/SelenzymeRF_scripts/filtered_dbs/MetaNetX_seed0

    python generate_selenzyme_db.py \
        --base_db_dir /path/to/unzipped/data_2023 \
        --test_tsv data/Splits-Rhea/Scaffold/seed_splits/seed0/test.tsv \
        --id_scheme rhea \
        --output_dir methods/SelenzymeRF/SelenzymeRF_scripts/filtered_dbs/Rhea_Scaffold_seed0
"""

import argparse
import os
import re
import shutil
import sys

import pandas as pd

# Reaction-keyed files that must be filtered to remove test-set reactions.
# (filename, separator, "first column is the MNXR reaction id" -- all three
# of them, per the docstring above.)
FILTERED_FILES = [
    ("reac_prop.tsv", "\t"),
    ("reac_seqs.tsv", "\t"),
    ("reac_smi.csv", ","),
]

# Files that are NOT reaction-specific (or, in the case of the two .npz
# fingerprint files, are reaction-indexed but only ever looked up for
# reaction ids that survive the reac_prop.tsv filter -- see module docstring).
# These are passed through unchanged (symlinked by default).
PASSTHROUGH_FILES = [
    "reac_xref.tsv",
    "seq_org.tsv",
    "org_lineage.csv",
    "FP_Morg.npz",
    "FP_MorgRF.npz",
    "rxn_consensus_20160612.txt",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a SelenzymeRF reference DB with a seed's test-set reactions excluded."
    )
    parser.add_argument("--base_db_dir", required=True,
                        help="Path to the unzipped data_2023 directory "
                             "(from compressed_data/data_2023.zip), containing "
                             "reac_prop.tsv, reac_seqs.tsv, reac_smi.csv, reac_xref.tsv, "
                             "seq_org.tsv, org_lineage.csv, FP_Morg.npz, FP_MorgRF.npz, "
                             "rxn_consensus_20160612.txt.")
    parser.add_argument("--test_tsv", required=True,
                        help="Path to the seed's split TSV whose reactions drive filtering -- "
                             "test.tsv for --filter_mode exclude, train.tsv for --filter_mode "
                             "include_only (flag name kept as --test_tsv for both since most "
                             "callers use the exclude mode; the id-resolution logic itself "
                             "doesn't care which split the file represents). "
                             "For --id_scheme metanetx: needs a 'reaction_id' column "
                             "with MNXR-style ids (data/Splits-DBs/MetaNetX/seed_splits/seedN/test.tsv). "
                             "For --id_scheme rhea: needs a 'REACTION_ID' column with plain "
                             "Rhea integer ids (data/Splits-Rhea/{Scaffold,Stratified}/seed_splits/seedN/test.tsv).")
    parser.add_argument("--filter_mode", choices=["exclude", "include_only"], default="exclude",
                        help="exclude (default): keep everything in --base_db_dir except "
                             "reactions in --test_tsv (broad reference pool, prevents leakage "
                             "by removing only the exact test reactions). "
                             "include_only: keep ONLY reactions in --test_tsv (intended for "
                             "--test_tsv pointing at train.tsv -- builds a self-contained, "
                             "train-only reference DB, methodologically consistent with how "
                             "SIMMER's reference DB is built purely from the train split. "
                             "Reference-pool coverage will be much smaller and some reactions "
                             "may end up with zero or few enzyme sequences, since Rhea's train "
                             "split alone often has far fewer curated sequences per reaction "
                             "than the full external database.)")
    parser.add_argument("--id_scheme", required=True, choices=["metanetx", "rhea"],
                        help="metanetx: test_tsv's reaction_id column already holds MNXR ids, "
                             "matched directly against column 1 of reac_prop.tsv. "
                             "rhea: test_tsv's REACTION_ID column holds plain Rhea integer ids; "
                             "these are resolved to MNXR ids via the 'rheaR:<id>' xref token in "
                             "reac_prop.tsv's own xref column (column 3), UNIONED with the "
                             "'rhea:<id> -> MNXR id' mappings in the standalone reac_xref.tsv "
                             "(needed because a single MNXR reaction can have more than one Rhea "
                             "id -- Rhea reactions come in groups of up to 4 ids: master/LR/RL/"
                             "bidirectional -- and reac_prop.tsv's xref column only ever stores one "
                             "of them; see report for details). Using the union is deliberately "
                             "over-inclusive: excluding a few extra non-test rows is harmless, "
                             "missing an actual test reaction is not.")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to write the filtered, self-contained DB into "
                             "(ready to be mounted/copied by start_server*.sh).")
    parser.add_argument("--seqs_fasta", default=None,
                        help="Optional path to seqs.fasta (unzipped from compressed_data/seqs.zip). "
                             "If given, it is symlinked/copied into --output_dir as well so the "
                             "output directory is fully self-contained.")
    parser.add_argument("--link_mode", choices=["symlink", "copy"], default="symlink",
                        help="How to place unfiltered passthrough files (reac_xref.tsv, seq_org.tsv, "
                             "org_lineage.csv, FP_Morg.npz, FP_MorgRF.npz, "
                             "rxn_consensus_20160612.txt, seqs.fasta) into --output_dir. "
                             "symlink (default) avoids duplicating multi-hundred-MB files per seed; "
                             "use copy if the output directory needs to be relocated/mounted somewhere "
                             "symlinks won't resolve (e.g. copied into a container image).")
    return parser.parse_args()


RHEA_XREF_RE = re.compile(r"^rheaR:(\d+)$")


def get_metanetx_ids_to_exclude(test_tsv: str) -> set:
    """metanetx scheme: test_tsv's reaction_id column already holds MNXR ids."""
    test_df = pd.read_csv(test_tsv, sep="\t")
    col = "reaction_id" if "reaction_id" in test_df.columns else "REACTION_ID"
    if col not in test_df.columns:
        raise ValueError(
            f"{test_tsv}: expected a 'reaction_id' column for --id_scheme metanetx, "
            f"found columns {list(test_df.columns)}"
        )
    return set(test_df[col].astype(str))


def get_rhea_ids_to_exclude(test_tsv: str, base_db_dir: str) -> set:
    """rhea scheme: resolve plain Rhea integer ids in test_tsv to MNXR ids.

    Uses both:
      1. reac_prop.tsv's own xref column (col 3, single 'rheaR:<id>' token per row), and
      2. the standalone reac_xref.tsv's 'rhea:<id>\tMNXR_id' rows (many-to-one: a single
         MNXR reaction can be tied to several Rhea ids -- the master id plus its LR/RL/
         bidirectional variants -- only one of which reac_prop.tsv's xref column records).
    The union of both is used so we don't accidentally leave a test reaction in the DB
    just because the test split cites a Rhea id variant that reac_prop.tsv doesn't store.
    """
    test_df = pd.read_csv(test_tsv, sep="\t")
    col = "REACTION_ID" if "REACTION_ID" in test_df.columns else "reaction_id"
    if col not in test_df.columns:
        raise ValueError(
            f"{test_tsv}: expected a 'REACTION_ID' column for --id_scheme rhea, "
            f"found columns {list(test_df.columns)}"
        )
    rhea_ids = set(test_df[col].astype(str))

    excluded_mnxr = set()

    # Source 1: reac_prop.tsv's own xref column (single-valued -- verified no ';'-separated
    # multi-values are ever present in this column; each row has exactly one xref).
    reac_prop_path = os.path.join(base_db_dir, "reac_prop.tsv")
    with open(reac_prop_path) as handler:
        for line in handler:
            if line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 3:
                continue
            mnxr_id, xref = row[0], row[2]
            m = RHEA_XREF_RE.match(xref)
            if m and m.group(1) in rhea_ids:
                excluded_mnxr.add(mnxr_id)

    # Source 2: standalone reac_xref.tsv, prefix 'rhea:' (not 'rheaR:'), format:
    #   rhea:<id>\tMNXR_id\t<extra columns...>
    reac_xref_path = os.path.join(base_db_dir, "reac_xref.tsv")
    if os.path.exists(reac_xref_path):
        with open(reac_xref_path) as handler:
            for line in handler:
                if line.startswith("#"):
                    continue
                row = line.rstrip("\n").split("\t")
                if len(row) < 2:
                    continue
                xref, mnxr_id = row[0], row[1]
                if xref.startswith("rhea:") and xref[len("rhea:"):] in rhea_ids:
                    excluded_mnxr.add(mnxr_id)
    else:
        print(f"WARNING: {reac_xref_path} not found; relying solely on reac_prop.tsv's own "
              f"xref column for rhea id matching (may under-exclude reactions that only "
              f"cite a non-master Rhea id variant).", file=sys.stderr)

    return excluded_mnxr


def filter_file(base_db_dir: str, output_dir: str, filename: str, sep: str,
                 ids: set, filter_mode: str) -> tuple:
    in_path = os.path.join(base_db_dir, filename)
    out_path = os.path.join(output_dir, filename)
    keep = (lambda first_col: first_col in ids) if filter_mode == "include_only" \
        else (lambda first_col: first_col not in ids)

    if filename == "reac_seqs.tsv":
        # No header, ~2.6M rows, 150+MB -- filter line-by-line rather than loading
        # the whole thing into a DataFrame.
        n_total, n_kept = 0, 0
        with open(in_path) as fin, open(out_path, "w") as fout:
            for line in fin:
                n_total += 1
                first_col = line.split("\t", 1)[0]
                if keep(first_col):
                    fout.write(line)
                    n_kept += 1
        return n_total, n_kept

    df = pd.read_csv(in_path, sep=sep, header=0 if filename == "reac_smi.csv" else None,
                      dtype=str, keep_default_na=False)
    first_col = df.columns[0]
    n_total = len(df)
    mask = df[first_col].isin(ids) if filter_mode == "include_only" else ~df[first_col].isin(ids)
    filtered_df = df[mask]
    n_kept = len(filtered_df)
    filtered_df.to_csv(out_path, sep=sep, index=(filename == "reac_smi.csv" and False),
                        header=(filename == "reac_smi.csv"))
    return n_total, n_kept


def link_or_copy(src: str, dst: str, link_mode: str) -> None:
    if os.path.lexists(dst):
        os.remove(dst)
    if link_mode == "symlink":
        os.symlink(os.path.abspath(src), dst)
    else:
        shutil.copy2(src, dst)


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    if args.id_scheme == "metanetx":
        ids = get_metanetx_ids_to_exclude(args.test_tsv)
    else:
        ids = get_rhea_ids_to_exclude(args.test_tsv, args.base_db_dir)

    verb = "inclusion (train-only DB)" if args.filter_mode == "include_only" else "exclusion"
    print(f"[{args.id_scheme}] {len(ids)} reaction id(s) resolved for {verb} "
          f"from {args.test_tsv}")
    if len(ids) == 0:
        print(f"WARNING: zero reaction ids resolved -- the filtered DB will be "
              f"{'EMPTY' if args.filter_mode == 'include_only' else 'IDENTICAL to the base DB (nothing excluded)'}. "
              "This almost certainly indicates an --id_scheme / column-name / id-format mismatch; "
              "double-check before using this output for evaluation.", file=sys.stderr)

    # 1. Filter the reaction-keyed lookup files.
    for filename, sep in FILTERED_FILES:
        in_path = os.path.join(args.base_db_dir, filename)
        if not os.path.exists(in_path):
            print(f"WARNING: {in_path} not found, skipping.", file=sys.stderr)
            continue
        n_total, n_kept = filter_file(args.base_db_dir, args.output_dir, filename, sep, ids, args.filter_mode)
        print(f"  {filename}: kept {n_kept}/{n_total} rows")

    # 2. Pass through the files that don't need filtering (see module docstring for why the
    #    two .npz fingerprint files are safe to leave untouched).
    for filename in PASSTHROUGH_FILES:
        src = os.path.join(args.base_db_dir, filename)
        if not os.path.exists(src):
            print(f"WARNING: {src} not found, skipping passthrough.", file=sys.stderr)
            continue
        dst = os.path.join(args.output_dir, filename)
        link_or_copy(src, dst, args.link_mode)
        print(f"  {filename}: {args.link_mode}ed unchanged")

    # 3. seqs.fasta (large, lives outside base_db_dir since it comes from a separate zip).
    if args.seqs_fasta:
        if os.path.exists(args.seqs_fasta):
            dst = os.path.join(args.output_dir, "seqs.fasta")
            link_or_copy(args.seqs_fasta, dst, args.link_mode)
            print(f"  seqs.fasta: {args.link_mode}ed unchanged from {args.seqs_fasta}")
        else:
            print(f"WARNING: --seqs_fasta {args.seqs_fasta} not found, skipping.", file=sys.stderr)

    print(f"\nFiltered SelenzymeRF DB written to {args.output_dir}")


if __name__ == "__main__":
    main()
