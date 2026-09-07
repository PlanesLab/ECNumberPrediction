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
parser.add_argument(
    "--pretrained_model", type=str, default=None,
    help="Path to pretrained model to finetune (default: rxnfp's shipped bert_pretrained -- "
         "a genuine MLM-pretrained-only checkpoint, NOT model/trained_512 or bert_class_ec_final, "
         "both of which are already fine-tuned BertForSequenceClassification checkpoints)",
)
parser.add_argument("--train_data", type=str, required=True, help="Dataset to finetune on")
parser.add_argument("--output_dir", type=str, required=True, help="Where to save the fine-tuned model")
parser.add_argument("--seed", type=int, default=42)
args_cli = parser.parse_args()

if args_cli.pretrained_model is None:
    args_cli.pretrained_model = pkg_resources.resource_filename("rxnfp", "models/transformers/bert_pretrained")

random.seed(args_cli.seed)

from rxnfp.models import SmilesClassificationModel
logger = logging.getLogger(__name__)

# from dotenv import load_dotenv, find_dotenv
# load_dotenv(find_dotenv())

df = pd.read_csv(args_cli.train_data)
print(df[['rxn', 'ec_subsubclass_label']].head())
train_df = df.loc[df['split']=='train']
print(train_df[['rxn', 'ec_subsubclass_label']].head())
eval_df = df[['rxn', 'class_id']].loc[df['split']=='val']
eval_df.columns = ['text', 'labels']
print(eval_df.head())

all_train_reactions = train_df.rxn.values.tolist()
corresponding_labels = train_df.class_id.values.tolist()
final_train_df = pd.DataFrame({'text': all_train_reactions, 'labels': corresponding_labels})
final_train_df = final_train_df.sample(frac=1., random_state=args_cli.seed)

model_args = {
    'wandb_project': None, 'num_train_epochs': 48, 'overwrite_output_dir': True,
    'learning_rate': 1e-5, 'gradient_accumulation_steps': 1,
    'regression': False, "num_labels": 353, "fp16": False,
    "evaluate_during_training": True, 'manual_seed': args_cli.seed,
    "max_seq_length": 512, "train_batch_size": 8,"warmup_ratio": 0.00,
    'output_dir': args_cli.output_dir,
    'thread_count': 4,
    # SmilesClassificationModel.__init__ first loads model_args.json from the --pretrained_model
    # checkpoint dir (self.args = self._load_model_args(model_name)), THEN applies this dict on
    # top via update_from_dict -- so this "labels_list": [] here overrides whatever labels_list
    # that checkpoint happened to save (e.g. model/trained_512's is 308 entries long, for an
    # unrelated external dataset/label space). Without this override, `if self.args.labels_list:
    # assert num_labels == len(self.args.labels_list)` hard-fails whenever num_labels=353 doesn't
    # match a warm-start checkpoint's own label count. Clearing it here makes the classifier head
    # size purely a function of num_labels, reinitialized fresh regardless of the checkpoint.
    # Also clear labels_map: leaving it unset here means it inherits trained_512's own
    # (stale, from an unrelated 308-class problem) mapping, since the labels_list==[] above
    # short-circuits the constructor code path that would otherwise rebuild/normalize it.
    "labels_list": [], "labels_map": None,
    }

model_path = args_cli.pretrained_model
print(model_path)
# ignore_mismatched_sizes threads through SmilesClassificationModel's **kwargs to the
# underlying from_pretrained() call -- without it, HF hard-errors on the classifier
# head's shape mismatch (308 saved vs 353 requested) instead of reinitializing it fresh.
model = SmilesClassificationModel("bert", model_path, num_labels=353, args=model_args, use_cuda=torch.cuda.is_available(), ignore_mismatched_sizes=True)


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