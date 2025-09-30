from pathlib import Path
import typer
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit
from drfp import DrfpEncoder

def main(path: str, output_path: str):
    # Load the data
    df = pd.read_csv(path)

    # Filter to rows where 'ec' is not missing
    df = df[df.ec.notna()]

    # Split the EC column into components, and create hierarchical labels
    df[["ec_1", "ec_2", "ec_3", "ec_4"]] = df.ec.str.split(".", expand=True)
    df["ec1"] = df.ec_1.astype(str)
    df["ec12"] = df.ec_1.astype(str) + "." + df.ec_2.astype(str)
    df["ec123"] = df.ec_1.astype(str) + "." + df.ec_2.astype(str) + "." + df.ec_3.astype(str)

    # Exclude reactions with first EC level equal to "7"
    df = df[df.ec1 != "7"]

    # Generate fingerprints using the DrfpEncoder
    df["fps"] = DrfpEncoder.encode(
        df.rxn,
        show_progress_bar=True,
        root_central_atom=False,
        radius=2,
        include_hydrogens=True,
        n_folded_length=10240,
    )

    # For each EC level (ec1, ec12, ec123) perform a stratified split for train, validation, and test
    for ec in ["ec1", "ec12", "ec123"]:
        X = df.rxn.to_numpy()
        y = df[ec].to_numpy()
        fps_arr = df.fps.to_numpy()
        groups = df.ec_1.to_numpy()  # stratify based on the first EC level

        sss_outer = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
        for i, (trainval_index, test_index) in enumerate(sss_outer.split(X, groups)):
            # Split train+val and test
            X_trainval, X_test = X[trainval_index], X[test_index]
            y_trainval, y_test = y[trainval_index], y[test_index]
            fps_trainval, fps_test = fps_arr[trainval_index], fps_arr[test_index]
            groups_trainval = groups[trainval_index]

            # Now split trainval into train and val
            sss_inner = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
            train_index, val_index = next(sss_inner.split(X_trainval, groups_trainval))

            X_train, X_val = X_trainval[train_index], X_trainval[val_index]
            y_train, y_val = y_trainval[train_index], y_trainval[val_index]
            fps_train, fps_val = fps_trainval[train_index], fps_trainval[val_index]

            # Save splits
            df_train = pd.DataFrame({
                "rxn_smiles": X_train,
                "label": y_train,
                "fps": [";".join(map(str, fp)) for fp in fps_train]
            })
            df_val = pd.DataFrame({
                "rxn_smiles": X_val,
                "label": y_val,
                "fps": [";".join(map(str, fp)) for fp in fps_val]
            })
            df_test = pd.DataFrame({
                "rxn_smiles": X_test,
                "label": y_test,
                "fps": [";".join(map(str, fp)) for fp in fps_test]
            })

            df_train.to_csv(f"{output_path}-{i}-{ec}-train.csv", index=False)
            df_val.to_csv(f"{output_path}-{i}-{ec}-val.csv", index=False)
            df_test.to_csv(f"{output_path}-{i}-{ec}-test.csv", index=False)

if __name__ == "__main__":
    typer.run(main)
