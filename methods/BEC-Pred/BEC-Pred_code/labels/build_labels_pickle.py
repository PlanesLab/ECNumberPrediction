"""
Build a {numeric_label: EC_string} pickle (the format label_assigner.py
expects) from generate_becpred_db.py's ec_class_labels.csv output.

Needed per-seed: the numeric label IDs generate_becpred_db.py assigns
depend on first-seen order over that seed's own train/test data, so a
fixed/stale labels_metanetx.pkl from a different run's data would silently
mismap predictions to the wrong EC numbers.
"""

import argparse
import pickle

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ec_class_labels_csv", required=True, help="ec_class_labels.csv from generate_becpred_db.py")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.ec_class_labels_csv)
    labels_dict = dict(zip(df["Assigned Label"], df["EC Class"]))

    with open(args.output, "wb") as f:
        pickle.dump(labels_dict, f)
    print(f"Wrote {len(labels_dict)} label mappings to '{args.output}'")


if __name__ == "__main__":
    main()
