#!/usr/bin/env python3
"""
============================================================
Majority EC Prediction Aggregator (Weighted Top-5)
============================================================

Computes Top-1 and Top-5 EC number predictions (Enzyme Commission numbers)
via majority voting across multiple prediction methods' outputs.

Features:
  ✅ Collapses ECs to 3rd hierarchical level (1.2.3.x → 1.2.3)
  ✅ Simple majority Top-1
  ✅ Weighted Top-5 (5→1 weight by rank)
  ✅ Works for one or all entities
  ✅ Choose specific or all methods
  ✅ Optional output CSV
"""

import argparse
import logging
import pandas as pd
from collections import Counter
from typing import Dict, List, Optional, Any


# ---------------------------------------------------------------------
# 🧩 Utility Functions
# ---------------------------------------------------------------------

def collapse_to_third_level(predictions: Any) -> Optional[str]:
    """Collapse EC predictions to the 3rd hierarchical level."""
    if pd.isna(predictions) or str(predictions).strip().lower() in {"", "nan"}:
        return None

    preds = str(predictions).split(";")
    collapsed_groups = []
    for group in preds:
        sub_preds = [x.strip() for x in group.split("|") if x.strip()]
        collapsed = {".".join(x.split(".")[:3]) for x in sub_preds}
        collapsed_groups.append("|".join(sorted(collapsed)))
    return ";".join(collapsed_groups)


def extract_top1(prediction: Optional[str]) -> List[str]:
    """Extract the top-1 EC predictions (first semicolon group)."""
    if not prediction:
        return []
    first_group = prediction.split(";")[0]
    return [x.strip() for x in first_group.split("|") if x.strip()]


def extract_top5(prediction: Optional[str]) -> List[str]:
    """Extract up to 5 EC predictions (flattening groups)."""
    if not prediction:
        return []
    groups = prediction.split(";")
    top5 = []
    for grp in groups:
        top5.extend(x.strip() for x in grp.split("|") if x.strip())
        if len(top5) >= 5:
            break
    return top5[:5]


def majority_vote_top1(preds_dict: Dict[str, List[str]]) -> Optional[str]:
    """Simple majority vote for Top-1 predictions."""
    all_preds = [p for preds in preds_dict.values() for p in preds]
    if not all_preds:
        return None
    return Counter(all_preds).most_common(1)[0][0]


def weighted_majority_vote_top5(preds_dict: Dict[str, List[str]]) -> Optional[str]:
    """
    Weighted majority vote for Top-5 predictions.
    Higher-ranked predictions get more weight (5→1).
    """
    weighted = Counter()
    for preds in preds_dict.values():
        for rank, p in enumerate(preds[:5]):  # limit to 5
            weight = max(5 - rank, 1)
            weighted[p] += weight

    if not weighted:
        return None
    return weighted.most_common(1)[0][0]


def compute_majority_for_row(row: pd.Series, methods: List[str]) -> pd.Series:
    """Compute Top-1 and weighted Top-5 majority votes for one row."""
    top1_dict = {m: extract_top1(row.get(m)) for m in methods}
    top5_dict = {m: extract_top5(row.get(m)) for m in methods}

    return pd.Series({
        "majority_vote_top1": majority_vote_top1(top1_dict),
        "majority_vote_top5": weighted_majority_vote_top5(top5_dict)
    })


# ---------------------------------------------------------------------
# 🚀 Main Logic
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute Top-1 and Weighted Top-5 majority EC predictions from multiple methods."
    )
    parser.add_argument("--input_csv", required=True, help="Path to merged prediction CSV.")
    parser.add_argument("--entity", required=True, help="Entity/reaction identifier (or 'all').")
    parser.add_argument("--output_csv", default=None, help="Optional output CSV path.")
    parser.add_argument("--methods", nargs="*", default=None, help="Specific method columns to use.")
    parser.add_argument("--use_all", action="store_true", help="Use all EC prediction columns automatically.")
    parser.add_argument("--id_col", default="entity", help="Name of the identifier column (default: 'entity').")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    # --- Load CSV ---
    try:
        df = pd.read_csv(args.input_csv)
    except FileNotFoundError:
        logging.error(f"File not found: {args.input_csv}")
        return

    df = df.fillna("nan")

    # --- Detect EC prediction columns ---
    excluded_cols = {args.id_col.lower(), "reaction", "rxn"}
    all_method_cols = [c for c in df.columns if c.lower() not in excluded_cols]

    # --- Determine which methods to use ---
    if args.use_all:
        methods = all_method_cols
        logging.info(f"Using ALL available EC prediction columns: {methods}")
    elif args.methods:
        methods = [m for m in args.methods if m in all_method_cols]
        if not methods:
            logging.error("No valid methods found in your selection.")
            return
        logging.info(f"Using SELECTED methods: {methods}")
    else:
        logging.error("No methods provided. Use --methods or --use_all.")
        return

    # --- Normalize EC formats ---
    for m in methods:
        df[m] = df[m].apply(collapse_to_third_level)

    # --- Compute results ---
    entity_col = args.id_col
    if entity_col not in df.columns:
        logging.error(f"Identifier column '{entity_col}' not found in data.")
        return

    if args.entity.lower() != "all":
        row_match = df[df[entity_col].astype(str).str.lower() == args.entity.lower()]
        if row_match.empty:
            logging.warning(f"No data found for entity: {args.entity}")
            return
        result = compute_majority_for_row(row_match.iloc[0], methods)
        print(f"\nEntity: {args.entity}")
        print(f"Top-1 Majority EC: {result['majority_vote_top1']}")
        print(f"Weighted Top-5 Majority EC: {result['majority_vote_top5']}")
    else:
        logging.info("Computing majority votes for ALL entities...")
        results = df.copy()
        results[["majority_vote_top1", "majority_vote_top5"]] = results.apply(
            compute_majority_for_row, axis=1, methods=methods
        )
        print(results[[entity_col, "majority_vote_top1", "majority_vote_top5"]].to_string(index=False))
        if args.output_csv:
            results[[entity_col, "majority_vote_top1", "majority_vote_top5"]].to_csv(args.output_csv, index=False)
            logging.info(f"Saved results to {args.output_csv}")


# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()
