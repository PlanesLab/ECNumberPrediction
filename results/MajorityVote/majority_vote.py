#!/usr/bin/env python3
"""
============================================================
Majority EC Prediction Aggregator (Weighted Top-N)
============================================================

Computes a ranked Top-N list of EC number predictions via weighted
majority voting across multiple prediction methods.

─────────────────────────────────────────────────────────────
SCORING  (example: --top_n 3)
─────────────────────────────────────────────────────────────
Each method's predictions are stored as groups (pipe-separated ECs
within a semicolon-delimited rank slot). All ECs in the same group
receive the same rank and the same points:

    rank-1 group  →  3 pts  per EC in that group
    rank-2 group  →  2 pts
    rank-3 group  →  1 pt

So if a method predicts "1.1.1|2.2.2" at rank-1, both 1.1.1 and
2.2.2 get 3 pts — they are not penalised for being co-predictions.

Points are summed across all methods → total weighted score per EC.

─────────────────────────────────────────────────────────────
RANKING & TIE-BREAKING  (all descending)
─────────────────────────────────────────────────────────────
    1. Total weighted score
    2. If tied → number of methods where this EC is in rank-1 group
    3. If still tied → number of methods where it is in rank-2 group
    4. … cascade through rank-N
    5. If still tied after all ranks → pipe-join into one shared slot
       e.g. "1.2.3|4.5.6" at position 3

─────────────────────────────────────────────────────────────
OUTPUT FORMAT
─────────────────────────────────────────────────────────────
  majority_top1          → best EC  (or "A|B" if unbreakable tie)
  majority_topN_ranked   → full ranked list, slots separated by " ; "
                           e.g. "1.1.1 ; 2.2.2 ; 3.3.3|4.4.4 ; 5.5.5"

─────────────────────────────────────────────────────────────
OTHER FEATURES
─────────────────────────────────────────────────────────────
  ✅ Collapses ECs to 3rd hierarchical level  (1.2.3.4 → 1.2.3)
  ✅ --top_n controls both scoring window and output length
  ✅ Groups (pipe-separated) at the same rank score equally
  ✅ Single entity or all entities
  ✅ Choose specific methods or --use_all
  ✅ Optional output CSV
"""

import argparse
import logging
import pandas as pd
from collections import Counter
from typing import Dict, List, Optional, Any, Tuple


# ---------------------------------------------------------------------
# 🧩 EC string parsing
# ---------------------------------------------------------------------

def collapse_to_third_level(predictions: Any) -> Optional[str]:
    """Collapse EC predictions to the 3rd hierarchical level (1.2.3.x → 1.2.3)."""
    if pd.isna(predictions) or str(predictions).strip().lower() in {"", "nan"}:
        return None
    collapsed_groups = []
    for group in str(predictions).split(";"):
        sub_preds = [x.strip() for x in group.split("|") if x.strip()]
        collapsed = {".".join(x.split(".")[:3]) for x in sub_preds}
        collapsed_groups.append("|".join(sorted(collapsed)))
    return ";".join(collapsed_groups)


def extract_topn_groups(prediction: Optional[str], n: int) -> List[List[str]]:
    """
    Parse a method's prediction string into an ordered list of rank-groups.

    Input format:  "1.1.1|2.2.2;3.3.3;4.4.4"
      - ';' separates rank slots
      - '|' separates co-predictions within the same rank slot

    Returns a list of groups, e.g.:
        [["1.1.1", "2.2.2"], ["3.3.3"], ["4.4.4"]]

    Stops accumulating once n ECs have been collected, but never splits
    a group mid-way — all members of a group are included or none are,
    since they share the same rank and cannot be ordered against each other.
    """
    if not prediction:
        return []
    groups = []
    total = 0
    for grp in prediction.split(";"):
        members = [x.strip() for x in grp.split("|") if x.strip()]
        if not members:
            continue
        groups.append(members)
        total += len(members)
        if total >= n:
            break
    return groups


# ---------------------------------------------------------------------
# 🧩 Scoring
# ---------------------------------------------------------------------

def build_score_tables(
    preds_dict: Dict[str, List[List[str]]], n: int
) -> Tuple[Counter, Dict[str, List[int]]]:
    """
    Compute score tables from all methods' grouped predictions.

    weighted_score[ec]:
        Sum of (n - rank) across all appearances, where rank is the
        0-indexed group position.  All ECs in the same group get the
        same points (n - rank), so co-predictions are treated fairly.

    topk_counts[ec]:
        List of length n.  topk_counts[ec][k] = number of methods that
        included this EC in their rank-k group.
        Used for cascade tie-breaking when weighted scores are equal.
    """
    weighted_score: Counter = Counter()
    topk_counts: Dict[str, List[int]] = {}

    for groups in preds_dict.values():
        for rank, members in enumerate(groups):
            pts = n - rank                      # same pts for every EC in this group
            for ec in members:
                weighted_score[ec] += pts
                if ec not in topk_counts:
                    topk_counts[ec] = [0] * n
                if rank < n:
                    topk_counts[ec][rank] += 1

    return weighted_score, topk_counts


# ---------------------------------------------------------------------
# 🧩 Ranking
# ---------------------------------------------------------------------

def rank_ecs(
    weighted_score: Counter,
    topk_counts: Dict[str, List[int]],
    n: int
) -> List[str]:
    """
    Produce a ranked list of up to n rank-slots.

    Sort key per EC (all descending):
        (weighted_score, count_at_rank0, count_at_rank1, ..., count_at_rank_(n-1))

    ECs identical on every dimension of this key are true ties and are
    pipe-joined into one shared slot.  Each slot counts as 1 toward n.
    """
    if not weighted_score:
        return []

    def sort_key(ec: str) -> tuple:
        counts = topk_counts.get(ec, [0] * n)
        return (weighted_score[ec], *counts)

    all_ecs = sorted(weighted_score.keys(), key=sort_key, reverse=True)

    ranked: List[str] = []
    i = 0
    while i < len(all_ecs) and len(ranked) < n:
        ec = all_ecs[i]
        key = sort_key(ec)
        tied = [e for e in all_ecs[i:] if sort_key(e) == key]
        ranked.append("|".join(sorted(tied)))
        i += len(tied)

    return ranked[:n]


# ---------------------------------------------------------------------
# 🧩 Per-row computation
# ---------------------------------------------------------------------

def compute_majority_for_row(row: pd.Series, methods: List[str], n: int) -> pd.Series:
    """Compute the majority-voted ranked Top-N list for one entity row."""
    preds_dict = {m: extract_topn_groups(row.get(m), n) for m in methods}
    weighted_score, topk_counts = build_score_tables(preds_dict, n)
    ranked = rank_ecs(weighted_score, topk_counts, n)

    top1     = ranked[0] if ranked else None
    topn_str = " ; ".join(ranked) if ranked else None

    return pd.Series({
        "majority_top1": top1,
        f"majority_top{n}_ranked": topn_str,
    })


# ---------------------------------------------------------------------
# 🚀 Main
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute a ranked Top-N majority EC list.\n\n"
            "Scoring: rank-1 = top_n pts, rank-2 = top_n-1 pts, ..., rank-N = 1 pt.\n"
            "Co-predictions (pipe-separated) at the same rank all receive equal points.\n"
            "Tie-breaking: weighted score → rank-1 count → rank-2 count → ... → pipe-join."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input_csv",  required=True, help="Path to merged prediction CSV.")
    parser.add_argument("--entity",     required=True, help="Entity identifier to query, or 'all'.")
    parser.add_argument("--output_csv", default=None,  help="Optional path to save results CSV.")
    parser.add_argument("--methods",    nargs="*",     help="Specific method columns to use.")
    parser.add_argument("--use_all",    action="store_true", help="Use all EC prediction columns.")
    parser.add_argument("--id_col",     default="entity", help="Name of the ID column (default: 'entity').")
    parser.add_argument(
        "--top_n", type=int, default=5,
        help=(
            "Number of ranked predictions to return (default: 5). "
            "Also sets the scoring window: rank-1 = top_n pts, ..., rank-N = 1 pt."
        )
    )
    args = parser.parse_args()

    if args.top_n < 1:
        parser.error("--top_n must be at least 1.")

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    try:
        df = pd.read_csv(args.input_csv)
    except FileNotFoundError:
        logging.error(f"File not found: {args.input_csv}")
        return

    df = df.fillna("nan")

    excluded = {args.id_col.lower(), "reaction", "rxn"}
    all_method_cols = [c for c in df.columns if c.lower() not in excluded]

    if args.use_all:
        methods = all_method_cols
        logging.info(f"Using ALL prediction columns: {methods}")
    elif args.methods:
        methods = [m for m in args.methods if m in all_method_cols]
        if not methods:
            logging.error("None of the specified methods were found in the CSV.")
            return
        logging.info(f"Using methods: {methods}")
    else:
        logging.error("Specify --methods or --use_all.")
        return

    for m in methods:
        df[m] = df[m].apply(collapse_to_third_level)

    entity_col = args.id_col
    if entity_col not in df.columns:
        logging.error(f"ID column '{entity_col}' not found.")
        return

    n = args.top_n
    topn_col = f"majority_top{n}_ranked"

    if args.entity.lower() != "all":
        row_match = df[df[entity_col].astype(str).str.lower() == args.entity.lower()]
        if row_match.empty:
            logging.warning(f"No data found for entity: {args.entity}")
            return

        result = compute_majority_for_row(row_match.iloc[0], methods, n)

        print(f"\nEntity : {args.entity}")
        print(f"top_n  : {n}  (rank-1 = {n} pts, rank-2 = {n-1} pts, ..., rank-{n} = 1 pt)")
        print(f"Top-1  : {result['majority_top1']}")
        print(f"\nTop-{n} ranked list:")
        if result[topn_col]:
            for i, slot in enumerate(result[topn_col].split(" ; "), 1):
                note = "  ← top-1" if i == 1 else ""
                print(f"  {i}. {slot}{note}")
        else:
            print("  (no predictions)")

    else:
        logging.info(f"Computing majority votes (top_n={n}) for ALL entities...")
        results = df.copy()
        results[["majority_top1", topn_col]] = results.apply(
            lambda row: compute_majority_for_row(row, methods, n), axis=1
        )
        print(results[[entity_col, "majority_top1", topn_col]].to_string(index=False))
        if args.output_csv:
            results[[entity_col, "majority_top1", topn_col]].to_csv(args.output_csv, index=False)
            logging.info(f"Saved results to {args.output_csv}")


# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()