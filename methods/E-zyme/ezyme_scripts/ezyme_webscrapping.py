import os
import requests
from bs4 import BeautifulSoup
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Run E-zyme for substrate-product pairs from a CSV")
    parser.add_argument("-i", "--input_csv", required=True, help="Path to input CSV file")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save results")
    parser.add_argument("--reaction_id_col", required=True, help="Column name for reaction ID")
    parser.add_argument("--reactant_col", required=True, help="Column name for reactant")
    parser.add_argument("--product_col", required=True, help="Column name for product")
    parser.add_argument("--delimiter", default=",", help="CSV delimiter (default ',')")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    reactions = pd.read_csv(args.input_csv, delimiter=args.delimiter)

    url = "https://www.genome.jp/tools-bin/e-zyme-ko"

    for _, row in reactions.iterrows():
        reaction_id = str(row[args.reaction_id_col])
        reactant = str(row[args.reactant_col])
        product = str(row[args.product_col])

        folder_path = os.path.join(args.output_dir, reaction_id)
        os.makedirs(folder_path, exist_ok=True)

        ezyme2_file = os.path.join(folder_path, f"{reactant}_{product}_E-zyme2.csv")
        ezyme1_file = os.path.join(folder_path, f"{reactant}_{product}_E-zyme1.csv")

        if os.path.exists(ezyme2_file) and os.path.exists(ezyme1_file):
            print(f"⏩ Skipping {reaction_id} – already processed.")
            continue

        # First POST request (mode=view)
        response1 = requests.post(url, data={"mode": "view", "cpd1": reactant, "cpd2": product})
        if response1.status_code != 200:
            print(f"❌ Failed to fetch data for {reaction_id} (Reactant: {reactant}, Product: {product})")
            continue

        soup1 = BeautifulSoup(response1.text, 'html.parser')
        try:
            id_input = soup1.find("input", {"name": "id"}).get("value")
        except Exception as e:
            print(f"❌ Couldn't extract ID for {reaction_id}: {e}")
            continue

        # Second POST request (mode=compute)
        response2 = requests.post(url, data={"mode": "compute", "cpd1": reactant, "cpd2": product, "id": id_input})
        soup2 = BeautifulSoup(response2.text, 'html.parser')

        try:
            # E-zyme2 table
            table1_div = soup2.find("div", id="ref_rp_img")
            if table1_div and table1_div.find("table"):
                table1 = table1_div.find("table")
                rows1 = table1.find_all("tr")
                data1 = []
                for r in rows1[1:]:
                    cols = r.find_all("td")
                    row_data = [col.get_text(separator=" ").strip().replace("\n", " ") for col in cols]
                    data1.append(row_data)
                df1 = pd.DataFrame(data1, columns=["RPAIR", "Score", "EC", "KO"])
                df1.to_csv(ezyme2_file, index=False)
                print(f"✅ Saved E-zyme2 table for {reaction_id}")
            else:
                print(f"⚠️ No E-zyme2 table found for {reaction_id}")

            # E-zyme1 table
            table2_div = soup2.find("div", id="ez1")
            if table2_div and table2_div.find("table"):
                table2 = table2_div.find("table")
                rows2 = table2.find_all("tr")
                data2 = []
                for r in rows2[1:]:
                    cols = r.find_all("td")
                    row_data = [col.get_text(separator=" ").strip().replace("\n", " ") for col in cols]
                    data2.append(row_data)
                df2 = pd.DataFrame(data2, columns=["EC", "Weighted Score", "Observed Freq", "Reactions"])
                df2.to_csv(ezyme1_file, index=False)
                print(f"✅ Saved E-zyme1 table for {reaction_id}")
            else:
                print(f"⚠️ No E-zyme1 table found for {reaction_id}")

        except Exception as e:
            print(f"❌ Error processing tables for {reaction_id}: {e}")

    print("✅ Processing complete.")

if __name__ == "__main__":
    main()
