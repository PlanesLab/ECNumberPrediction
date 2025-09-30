"""
Script: query_claire.py
Author: Josefina Arcagni
Date: 2025-09-11

Description:
    This script performs EC Number prediction for chemical reactions using the CLAIRE method.
    It loads test and training data, reaction IDs, and model files, then runs inference to predict EC numbers.
    The script is designed to be executed from the command line with specified input file paths and parameters.

Args:
    --test_data_path (str): Path to the test data (.npy file) containing reaction fingerprints.
    --train_data_path (str): Path to the training data (.pkl file) containing model lookup fingerprints.
    --train_labels_path (str): Path to the training labels (.pkl file) containing EC numbers.
    --test_csv_path (str): Path to the CSV file containing reaction IDs for the test set.
    --reaction_id_col (str): Column name in the CSV file that contains reaction IDs.
    --model_path (str): Path to the pretrained model (.pth file).
    --gmm_path (str): Path to the GMM ensemble (.pkl file).

Returns:
    Saves prediction results for the test set, including top-k EC number predictions per reaction.
"""
import pickle
import numpy as np
import pandas as pd
from dev.prediction.inference_EC import inference
import sys
import argparse
parser = argparse.ArgumentParser(description="CLAIRE EC Number Prediction")

parser.add_argument('--test_data_path', type=str, required=True, help='Path to test_data.npy')
parser.add_argument('--train_data_path', type=str, required=True, help='Path to model_lookup_train.pkl')
parser.add_argument('--train_labels_path', type=str, required=True, help='Path to labels_train_ec3.pkl')
parser.add_argument('--test_csv_path', type=str, required=True, help='Path to test CSV file')
parser.add_argument('--reaction_id_col', type=str, required=True, help='Reaction ID column name in CSV')
parser.add_argument('--model_path', type=str, required=True, help='Path to model (.pth)')
parser.add_argument('--gmm_path', type=str, required=True, help='Path to GMM ensemble (.pkl)')

args = parser.parse_args()

test_data_path = args.test_data_path
train_data_path = args.train_data_path
train_labels_path = args.train_labels_path
test_csv_path = args.test_csv_path
reaction_id_col = args.reaction_id_col
pretrained_model_path = args.pretrained_model_path
gmm_path = args.gmm_path

# Load the concatenated fingerprints (test data) from the .npy file
test_data = np.load(test_data_path)

# Load training data and labels
train_data = pickle.load(open(train_data_path, 'rb'))
train_labels = pickle.load(open(train_labels_path, 'rb'))

test_labels = None

# Load the Reaction IDs from the CSV file and use them as test tags
reaction_df = pd.read_csv(test_csv_path)
test_tags = reaction_df[reaction_id_col].tolist()

# Define the pretrained model and GMM ensemble paths
pretrained_model = pretrained_model_path

results = inference(train_data, test_data, train_labels, test_tags, test_labels, pretrained_model, evaluation=False, topk=3, gmm=gmm_path)

# results will appear in test_predictions 
