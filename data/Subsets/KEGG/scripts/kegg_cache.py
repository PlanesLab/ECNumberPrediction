"""
Shared reaction-level prediction cache for the KEGG bootstrap (Case 1, full-pool
sampling with replacement). The 10 seeds draw from the same 7780-reaction pool
and overlap heavily, so this lets a method skip any reaction it already
predicted in an earlier seed instead of re-querying it.

One cache file per method at results/Case1/kegg_pool_cache/<Method>.csv,
normalized to columns `Reaction ID,<pred_col>`. E-zyme is excluded -- it has
its own pair-level scrape cache (methods/E-zyme/output/outputKEGG_pool).

Subcommands:
  seed    Import an existing results file (original or a past seed's raw output)
          into a method's cache.
  split   Given a seed's full sampled reactions CSV + a method's cache, write the
          "todo" subset (rows not yet cached) as a drop-in replacement
          kegg_reactions_current_test.csv + 8:2KEGGTest_canonicalized.txt pair.
  merge   Given a seed's full sampled reactions CSV + a method's fresh raw output
          for the todo subset (+ the todo CSV, needed if the fresh output has no ID
          column of its own), append the fresh predictions to the cache and write
          the seed's full per-method result (cache lookups re-expand any duplicate
          draws from the with-replacement sample).
"""

import argparse
import os
import sys

import pandas as pd
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'Preprocessing'))
from canonicalize_rxn_SMILES import canonicalize_reaction_smiles  # noqa: E402

CACHE_DIR = "results/Case1/kegg_pool_cache"


def cache_path(method: str) -> str:
    return os.path.join(CACHE_DIR, f"{method}.csv")


def load_cache(method: str) -> pd.DataFrame:
    path = cache_path(method)
    if os.path.exists(path):
        return pd.read_csv(path, dtype=str).fillna("")
    return pd.DataFrame(columns=["Reaction ID", "pred"])


def save_cache(method: str, df: pd.DataFrame) -> None:
    os.makedirs(CACHE_DIR, exist_ok=True)
    df.to_csv(cache_path(method), index=False)


def cmd_seed(args: argparse.Namespace) -> None:
    src = pd.read_csv(args.input, dtype=str).fillna("")
    if args.positional_ids_from:
        ids_df = pd.read_csv(args.positional_ids_from, dtype=str)
        if len(ids_df) != len(src):
            print(f"WARNING: {args.input} has {len(src)} rows but {args.positional_ids_from} "
                  f"has {len(ids_df)}; truncating to the shorter length")
            n = min(len(ids_df), len(src))
            src = src.iloc[:n].reset_index(drop=True)
            ids_df = ids_df.iloc[:n].reset_index(drop=True)
        rxn_ids = ids_df[args.id_col]
    else:
        rxn_ids = src[args.id_col]

    new_rows = pd.DataFrame({"Reaction ID": rxn_ids, "pred": src[args.pred_col]})
    new_rows = new_rows.drop_duplicates(subset="Reaction ID", keep="first")

    cache = load_cache(args.method)
    combined = pd.concat([cache, new_rows], ignore_index=True)
    combined = combined.drop_duplicates(subset="Reaction ID", keep="first")
    save_cache(args.method, combined)
    print(f"{args.method}: cache now has {len(combined)} reactions "
          f"(+{len(combined) - len(cache)} new from {args.input})")


def cmd_split(args: argparse.Namespace) -> None:
    seed_sample = pd.read_csv(args.seed_sample, dtype=str).fillna("")
    cache = load_cache(args.method)
    cached_ids = set(cache["Reaction ID"])

    todo_mask = ~seed_sample["Reaction ID"].isin(cached_ids)
    todo = seed_sample[todo_mask].drop_duplicates(subset="Reaction ID").reset_index(drop=True)

    n_total = seed_sample["Reaction ID"].nunique()
    print(f"{args.method}: {n_total} unique reactions in seed sample, "
          f"{n_total - len(todo)} already cached, {len(todo)} to query")

    os.makedirs(args.output_dir, exist_ok=True)
    csv_path = os.path.join(args.output_dir, "kegg_reactions_current_test.csv")
    txt_path = os.path.join(args.output_dir, "8:2KEGGTest_canonicalized.txt")
    todo.to_csv(csv_path, index=False)

    if len(todo) > 0:
        canon = todo["Isomeric_SMILES"].astype(str).apply(canonicalize_reaction_smiles)
        with open(txt_path, "w") as f:
            for s in canon.fillna(""):
                f.write(s + "\n")
    else:
        open(txt_path, "w").close()

    print(f"Wrote {len(todo)} rows to '{csv_path}' and '{txt_path}'")


def cmd_merge(args: argparse.Namespace) -> None:
    if os.path.getsize(args.fresh_output) == 0:
        fresh_ids, fresh_preds = [], []
    else:
        fresh = pd.read_csv(args.fresh_output, dtype=str).fillna("")
        if args.positional:
            todo = pd.read_csv(args.todo_csv, dtype=str)
            n = min(len(todo), len(fresh))
            fresh_ids = todo["Reaction ID"].iloc[:n].tolist()
            fresh_preds = fresh[args.pred_col].iloc[:n].tolist()
        else:
            fresh_ids = fresh[args.id_col].tolist()
            fresh_preds = fresh[args.pred_col].tolist()

    new_rows = pd.DataFrame({"Reaction ID": fresh_ids, "pred": fresh_preds})
    new_rows = new_rows.drop_duplicates(subset="Reaction ID", keep="first")

    cache = load_cache(args.method)
    combined = pd.concat([cache, new_rows], ignore_index=True)
    combined = combined.drop_duplicates(subset="Reaction ID", keep="first")
    save_cache(args.method, combined)
    print(f"{args.method}: added {len(new_rows)} fresh predictions, cache now has {len(combined)} reactions")

    seed_sample = pd.read_csv(args.seed_sample, dtype=str).fillna("")
    merged = seed_sample[["Reaction ID"]].merge(combined, on="Reaction ID", how="left")
    n_missing = merged["pred"].isna().sum()
    if n_missing > 0:
        print(f"WARNING: {n_missing}/{len(merged)} rows in the seed sample still have no "
              f"cached prediction after this merge (method may have failed on them)")
    merged = merged.rename(columns={"pred": args.out_pred_col})
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    merged.to_csv(args.output, index=False)
    print(f"Wrote seed's full {len(merged)}-row result to '{args.output}'")


def main() -> None:
    parser = argparse.ArgumentParser(description="KEGG bootstrap reaction-level prediction cache.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_seed = sub.add_parser("seed")
    p_seed.add_argument("--method", required=True)
    p_seed.add_argument("--input", required=True)
    p_seed.add_argument("--id_col", default="Reaction ID")
    p_seed.add_argument("--pred_col", required=True)
    p_seed.add_argument("--positional_ids_from", default=None,
                         help="If the input has no usable ID column, pull IDs positionally from this CSV instead.")
    p_seed.set_defaults(func=cmd_seed)

    p_split = sub.add_parser("split")
    p_split.add_argument("--method", required=True)
    p_split.add_argument("--seed_sample", required=True)
    p_split.add_argument("--output_dir", required=True)
    p_split.set_defaults(func=cmd_split)

    p_merge = sub.add_parser("merge")
    p_merge.add_argument("--method", required=True)
    p_merge.add_argument("--seed_sample", required=True)
    p_merge.add_argument("--fresh_output", required=True)
    p_merge.add_argument("--id_col", default="Reaction ID")
    p_merge.add_argument("--pred_col", required=True)
    p_merge.add_argument("--positional", action="store_true",
                          help="Fresh output has no ID column; pair it positionally with --todo_csv instead.")
    p_merge.add_argument("--todo_csv", default=None)
    p_merge.add_argument("--out_pred_col", required=True)
    p_merge.add_argument("--output", required=True)
    p_merge.set_defaults(func=cmd_merge)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
