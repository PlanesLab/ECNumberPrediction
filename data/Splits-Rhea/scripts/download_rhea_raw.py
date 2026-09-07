"""
Downloads Rhea's raw distribution files needed to rebuild master.tsv from
source.

Fetches from the Expasy mirror (https://ftp.expasy.org/databases/rhea/tsv/):
  - rhea-directions.tsv       RHEA_ID_MASTER, RHEA_ID_LR, RHEA_ID_RL, RHEA_ID_BI
  - rhea2ec.tsv                RHEA_ID, DIRECTION, MASTER_ID, ID (EC number)
  - rhea-reaction-smiles.tsv   RHEA_ID, SMILES (LR/RL ids only)

Writes them unmodified into data/Splits-Rhea/raw/.
"""

import argparse
import os
import urllib.request

BASE_URL = "https://ftp.expasy.org/databases/rhea/tsv/"
FILES = ["rhea-directions.tsv", "rhea2ec.tsv", "rhea-reaction-smiles.tsv"]
OUTPUT_DIR = "data/Splits-Rhea/raw"


def main() -> None:
    parser = argparse.ArgumentParser(description="Download raw Rhea TSV distribution files.")
    parser.add_argument("--output_dir", default=OUTPUT_DIR)
    parser.add_argument("--force", action="store_true", help="Re-download even if the file already exists.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    for fname in FILES:
        dest = os.path.join(args.output_dir, fname)
        if os.path.exists(dest) and not args.force:
            print(f"Skipping {fname} (already present, use --force to re-download)")
            continue
        url = BASE_URL + fname
        print(f"Downloading {url} -> {dest}")
        urllib.request.urlretrieve(url, dest)
    print("Done.")


if __name__ == "__main__":
    main()
