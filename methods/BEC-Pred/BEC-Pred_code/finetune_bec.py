import os
import numpy as np
import pandas as pd
import torch
import logging
import random
import pkg_resources
import sklearn
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--pretrained_model", type=str, required=True, help="Path to pretrained model to finetune")
parser.add_argument("--train_data", type=str, required=True, help="Dataset to finetune on")
parser.add_argument("--output_dir", type=str, required=True, help="Where to save the fine-tuned model")
args_cli = parser.parse_args()

from rxnfp.models import SmilesClassificationModel
logger = logging.getLogger(__name__)

# from dotenv import load_dotenv, find_dotenv
# load_dotenv(find_dotenv())

df = pd.read_csv(args.train_data)
print(df[['rxn', 'ec_subsubclass_label']].head())
train_df = df.loc[df['split']=='train']
print(train_df[['rxn', 'ec_subsubclass_label']].head())
eval_df = df[['rxn', 'ec_subsubclass_label']].loc[df['split']=='val']
eval_df.columns = ['text', 'labels']
print(eval_df.head())

all_train_reactions = train_df.rxn.values.tolist()
corresponding_labels = train_df.class_id.values.tolist()
final_train_df = pd.DataFrame({'text': all_train_reactions, 'labels': corresponding_labels})
final_train_df = final_train_df.sample(frac=1.)

model_args = {
    'wandb_project': None, 'num_train_epochs': 48, 'overwrite_output_dir': True,
    'learning_rate': 1e-5, 'gradient_accumulation_steps': 1,
    'regression': False, "num_labels": 353, "fp16": False,
    "evaluate_during_training": True, 'manual_seed': 42,
    "max_seq_length": 512, "train_batch_size": 8,"warmup_ratio": 0.00,
    'output_dir': args_cli.output_dir,
    'thread_count': 4,
    }

# optional
model_path =  args_cli.input_dir
print(model_path)
model = SmilesClassificationModel("bert", model_path, num_labels=353, args=model_args, use_cuda=torch.cuda.is_available())


# optional
# train_model_path =  pkg_resources.resource_filename("best_model")

def f1_multiclass(labels, preds):
      return sklearn.metrics.f1_score(labels, preds, average='weighted')

def prec_multiclass(labels, preds):
      return sklearn.metrics.precision_score(labels, preds, average='weighted')

def rec_multiclass(labels, preds):
      return sklearn.metrics.recall_score(labels, preds, average='weighted')

model.train_model(final_train_df, eval_df=eval_df, prec=prec_multiclass, rec=rec_multiclass, acc=sklearn.metrics.accuracy_score, mcc=sklearn.metrics.matthews_corrcoef, f1=f1_multiclass)

model_output_dir = model.args.output_dir
print(f"Model output will be saved to: {model_output_dir}")

result, model_outputs, wrong_predictions = model.eval_model(eval_df, prec=prec_multiclass, rec=rec_multiclass, acc=sklearn.metrics.accuracy_score, mcc=sklearn.metrics.matthews_corrcoef, f1=f1_multiclass)