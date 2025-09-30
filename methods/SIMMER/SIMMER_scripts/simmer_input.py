import pandas as pd
import typer

def process_reactions(
    input_path: str,
    output_path: str,
    reaction_id_col: str = "reaction_id",
    substrates_col: str = "substrates_products",
    smiles_col: str = "reaction_smiles",
    ec_col: str = "EC_number",   # <-- EC column name is now fully configurable
    sep: str = ",",              # Separator for input file (CSV/TSV)
    reaction_sep: str = ">>",    # Separator between substrates and products
    include_ec: bool = False,
):
    # Load the input file
    df = pd.read_csv(input_path, sep=sep)

    output_df = pd.DataFrame()
    output_df["reaction"] = df[reaction_id_col]

    # Substrates and products
    output_df["left_comp"] = df[substrates_col].apply(
        lambda x: x.split(reaction_sep)[0].strip().replace(".", "//").replace('"', "").strip()
    )
    output_df["right_comp"] = df[substrates_col].apply(
        lambda x: x.split(reaction_sep)[1].strip().replace(".", "//").replace('"', "").strip()
    )

    # Split SMILES
    output_df[["left_smiles", "right_smiles"]] = df[smiles_col].str.split(">>", expand=True)
    output_df["left_smiles"] = output_df["left_smiles"].str.strip()
    output_df["right_smiles"] = output_df["right_smiles"].str.strip()

    output_df = output_df[
        output_df[["left_comp", "right_comp", "left_smiles", "right_smiles"]]
        .apply(lambda row: all(row.str.strip() != ""), axis=1)
    ]

    if include_ec:
        if ec_col in df.columns:
            output_df[ec_col] = df[ec_col]
        else:
            raise ValueError(f"⚠️ EC column '{ec_col}' not found in input file.")

    output_df.to_csv(output_path, index=False)
    print(f"✅ Conversion complete. Saved to {output_path}")


def main():
    typer.run(process_reactions)


if __name__ == "__main__":
    main()
