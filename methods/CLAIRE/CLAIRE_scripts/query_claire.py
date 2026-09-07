"""
Script: query_claire.py
Author: Josefina Arcagni
Date: 2025-09-11

Predicts EC numbers for reactions via CLAIRE: loads test/train fingerprints,
labels, and model files, then runs inference and saves top-k predictions
per reaction. See --help for arguments.
"""
import pickle
import numpy as np
import pandas as pd
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'CLAIRE_code', 'CLAIRE'))
from dev.prediction.inference_EC import inference
import argparse
parser = argparse.ArgumentParser(description="CLAIRE EC Number Prediction")

parser.add_argument('--test_data_path', type=str, required=True, help='Path to test_data.npy')
parser.add_argument('--train_data_path', type=str, required=True, help='Path to model_lookup_train.pkl')
parser.add_argument('--train_labels_path', type=str, required=True, help='Path to labels_train_ec3.pkl')
parser.add_argument('--test_csv_path', type=str, required=True, help='Path to test CSV file')
parser.add_argument('--reaction_id_col', type=str, required=True, help='Reaction ID column name in CSV')
parser.add_argument('--model_path', type=str, required=True, help='Path to model (.pth)')
parser.add_argument('--gmm_path', type=str, required=True, help='Path to GMM ensemble (.pkl)')
parser.add_argument('--out_filename', type=str, required=True,
                     help="Output path prefix WITHOUT extension -- inference() writes '<out_filename>_prediction.csv'")

args = parser.parse_args()

test_data_path = args.test_data_path
train_data_path = args.train_data_path
train_labels_path = args.train_labels_path
test_csv_path = args.test_csv_path
reaction_id_col = args.reaction_id_col
pretrained_model_path = args.model_path
gmm_path = args.gmm_path

# Load the concatenated fingerprints (test data) from the .npy file
test_data = np.load(test_data_path)

# Load training data and labels
train_data = pickle.load(open(train_data_path, 'rb'))
train_labels = pickle.load(open(train_labels_path, 'rb'))

test_labels = None

# Load the Reaction IDs from the CSV file and use them as test tags
reaction_df = pd.read_csv(test_csv_path, sep=None, engine='python')
test_tags = reaction_df[reaction_id_col].tolist()

# Define the pretrained model and GMM ensemble paths
pretrained_model = pretrained_model_path

os.makedirs(os.path.dirname(os.path.abspath(args.out_filename)), exist_ok=True)
results = inference(train_data, test_data, train_labels, test_tags, test_labels, pretrained_model,
                     evaluation=False, out_filename=args.out_filename, topk=3, gmm=gmm_path)
print(f"Predictions written to '{args.out_filename}_prediction.csv'") 
