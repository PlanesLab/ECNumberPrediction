# Script: create_fps.py
# Author: Josefina Arcagni
# Date: 2025-09-11
# Description: This script creates combined fingerprints for reactions using RXNFP and DRFP.

import sys
import os
import pickle
import numpy as np
import argparse
from rxnfp.transformer_fingerprints import (
    RXNBERTFingerprintGenerator, get_default_model_and_tokenizer
)

# -------------------------------
# Parse arguments
# -------------------------------

parser = argparse.ArgumentParser(description="Create combined fingerprints for reactions.")
parser.add_argument("query_path", type=str, help="Path to query reactions file")
parser.add_argument("fps_path", type=str, help="Path to DRFP fingerprints pickle file")
parser.add_argument("test_out", type=str, help="Output path for combined fingerprints (.npy)")

args = parser.parse_args()

query_path = args.query_path
fps_path = args.fps_path
test_out = args.test_out

# -------------------------------
# Setup model + tokenizer
# -------------------------------
model, tokenizer = get_default_model_and_tokenizer()
rxnfp_generator = RXNBERTFingerprintGenerator(model, tokenizer)

# -------------------------------
# Load query reactions
# -------------------------------
with open(query_path, "r") as f:
    rxns = [line.strip() for line in f if line.strip()]

# -------------------------------
# Generate RXNFP embeddings
# -------------------------------
rxnfp_out = os.path.join(os.path.dirname(test_out), "rxnfp_emb.pkl")
rxnfp = rxnfp_generator.convert_batch(rxns)
with open(rxnfp_out, "wb") as f:
    pickle.dump(rxnfp, f)

# -------------------------------
# Load DRFP + RXNFP fingerprints
# -------------------------------
with open(fps_path, "rb") as f:
    drfp = pickle.load(f)
with open(rxnfp_out, "rb") as f:
    rxnfp = pickle.load(f)

# -------------------------------
# Concatenate fingerprints
# -------------------------------
test_data = np.concatenate([
    np.concatenate((np.reshape(r, (1, 256)), np.reshape(d, (1, 256))), axis=1)
    for r, d in zip(rxnfp, drfp)
], axis=0)

# -------------------------------
# Save combined array
# -------------------------------
np.save(test_out, test_data)