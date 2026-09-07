"""
Derives a single substrate/product KEGG-compound-ID pair per reaction, for
E-zyme's webscraper (accepts only one substrate + one product per query).

Parses the `Equation` column, using `R Class`'s compound pair where available
to resolve which side is substrate vs product (R Class's pair members are
sorted alphabetically, not ordered by side).

Falls back to stripping common cofactors (water, ATP/ADP, NAD(P)(H), CoA,
Pi/PPi, O2, CO2, H+, ...) and taking the first remaining compound per side.
Reactions left empty on either side after filtering are dropped and reported.
"""

import argparse
import re

import pandas as pd

COFACTORS = {
    "C00001",  # water
    "C00002", "C00008", "C00020",  # ATP, ADP, AMP
    "C00003", "C00004",  # NAD+, NADH
    "C00005", "C00006",  # NADPH, NADP+
    "C00007",  # O2
    "C00009", "C00013",  # phosphate, pyrophosphate
    "C00010",  # CoA
    "C00011",  # CO2
    "C00080",  # H+
    "C00238", "C00305",  # K+, Mg2+
}


def parse_side(side: str):
    ids = []
    for term in side.split("+"):
        m = re.search(r"[A-Z]\d{5}", term)
        if m:
            ids.append(m.group(0))
    return ids


def pick_pair(equation: str, r_class):
    try:
        lhs, rhs = equation.split("<=>")
    except ValueError:
        return None, None
    lhs_ids, rhs_ids = parse_side(lhs), parse_side(rhs)

    if isinstance(r_class, str) and r_class.strip():
        m = re.search(r"([A-Z]\d{5})_([A-Z]\d{5})", r_class)
        if m:
            a, b = m.group(1), m.group(2)
            if a in lhs_ids and b in rhs_ids:
                return a, b
            if b in lhs_ids and a in rhs_ids:
                return b, a

    lhs_main = [c for c in lhs_ids if c not in COFACTORS]
    rhs_main = [c for c in rhs_ids if c not in COFACTORS]
    if lhs_main and rhs_main:
        return lhs_main[0], rhs_main[0]
    return None, None


def main() -> None:
    parser = argparse.ArgumentParser(description="Derive substrate/product compound-ID pairs for E-zyme.")
    parser.add_argument("--input", required=True, help="Reaction CSV (needs Reaction ID, Equation, R Class columns)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, dtype=str)
    reactants, products = [], []
    n_dropped = 0
    for _, row in df.iterrows():
        r_class = row["R Class"] if "R Class" in df.columns else None
        sub, prod = pick_pair(row["Equation"], r_class)
        if sub is None:
            n_dropped += 1
        reactants.append(sub)
        products.append(prod)

    df["Reactants"] = reactants
    df["Products"] = products
    out_df = df.dropna(subset=["Reactants", "Products"])[["Reaction ID", "Reactants", "Products"]]
    out_df.to_csv(args.output, index=False)
    print(f"Derived pairs for {len(out_df)}/{len(df)} reactions ({n_dropped} dropped, no usable substrate/product pair)")
    print(f"Wrote '{args.output}'")


if __name__ == "__main__":
    main()
