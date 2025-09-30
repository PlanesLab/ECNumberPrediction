# Script: get_ec_predictions.py
# Author: Josefina Arcagni
# Date: 2025-09-11
# Description: This script parses EC predictions from CLAIRE results.

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Parse EC predictions from CLAIRE results.")
    parser.add_argument("--input", required=True, help="Path to input CSV file")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    args = parser.parse_args()

    data = []
    with open(args.input, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            reaction_id = parts[0].strip()
            predictions = []
            for pred in parts[1:6]:
                pred_ec = pred.strip().split('/')[0]
                if pred_ec.startswith("EC:"):
                    pred_ec = pred_ec.replace("EC:", "", 1)
                predictions.append(pred_ec)
            all_predictions = ";".join(predictions)
            data.append({
                "reactionID": reaction_id,
                "predictions": all_predictions
            })

    df = pd.DataFrame(data)
    df.to_csv(args.output, index=False)
    print(f"Data saved to {args.output}")

if __name__ == "__main__":
    main()
