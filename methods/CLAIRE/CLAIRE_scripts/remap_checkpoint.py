"""
Remaps a train-pred_rxn_EC.py checkpoint's state_dict keys to the naming
inference_EC.py's LayerNormNet expects.

The two files define separate LayerNormNet classes for the same 5-layer MLP,
with different key naming ("layers.0.0.weight" vs "fc1.weight"), so a training
checkpoint doesn't load directly into the inference model.

Assumes --num_layers 5: layers.0 -> fc1/ln1, layers.1 -> fc2/ln2,
layers.2 -> fc4/ln4, layers.3 -> fc5/ln5, layers.4 -> fc3 (final plain Linear).
"""

import argparse

import torch

KEY_MAP = {
    "layers.0.0.weight": "fc1.weight", "layers.0.0.bias": "fc1.bias",
    "layers.0.1.weight": "ln1.weight", "layers.0.1.bias": "ln1.bias",
    "layers.1.0.weight": "fc2.weight", "layers.1.0.bias": "fc2.bias",
    "layers.1.1.weight": "ln2.weight", "layers.1.1.bias": "ln2.bias",
    "layers.2.0.weight": "fc4.weight", "layers.2.0.bias": "fc4.bias",
    "layers.2.1.weight": "ln4.weight", "layers.2.1.bias": "ln4.bias",
    "layers.3.0.weight": "fc5.weight", "layers.3.0.bias": "fc5.bias",
    "layers.3.1.weight": "ln5.weight", "layers.3.1.bias": "ln5.bias",
    "layers.4.weight": "fc3.weight", "layers.4.bias": "fc3.bias",
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Remap a CLAIRE training checkpoint for inference_EC.py's LayerNormNet.")
    parser.add_argument("--input", required=True, help="Checkpoint written by train-pred_rxn_EC.py (train_final.pth)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    state_dict = torch.load(args.input, map_location="cpu")
    missing = set(state_dict) - set(KEY_MAP)
    if missing:
        raise SystemExit(
            f"Unrecognized key(s) in {args.input}: {sorted(missing)} -- checkpoint doesn't match "
            "the expected num_layers=5 architecture this remapping assumes."
        )

    remapped = {KEY_MAP[k]: v for k, v in state_dict.items()}
    torch.save(remapped, args.output)
    print(f"Remapped {len(remapped)} params from '{args.input}' -> '{args.output}'")


if __name__ == "__main__":
    main()
