import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

def get_ec_subsubclass(ec_number):
    """Extract subsubclass from EC number (e.g., 3.5.1.50 -> 3.5.1)"""
    try:
        parts = str(ec_number).split('.')
        if len(parts) >= 3:
            return '.'.join(parts[:3])
        return str(ec_number)
    except:
        return None

# Load data
print("Loading data...")
df = pd.read_csv('/scratch/jarcagniriv/Case1-Rhea/Database/smiles_vs_ec_with_id.tsv', sep='\t')
print(f"Total reactions: {len(df)}")

# Extract EC subsubclass
print("\nExtracting EC subsubclasses...")
df['ec_subsubclass'] = df['EC_NUMBER'].apply(get_ec_subsubclass)
df_clean = df.dropna(subset=['ec_subsubclass'])
print(f"Reactions with valid EC: {len(df_clean)}")

# Stratified split by EC subsubclass
print("\nPerforming stratified random split (90% train, 10% val)...")
print("-" * 80)

train_df, val_df = train_test_split(
    df_clean,
    test_size=0.1,
    stratify=df_clean['ec_subsubclass'],
    random_state=42
)

# Statistics per EC subsubclass
ec_groups = df_clean.groupby('ec_subsubclass')
print("\nPer EC subsubclass distribution:")

for ec_class, group in ec_groups:
    n_total = len(group)
    n_train = len(train_df[train_df['ec_subsubclass'] == ec_class])
    n_val = len(val_df[val_df['ec_subsubclass'] == ec_class])
    
    print(f"EC {ec_class}: total {n_total:5d} → "
          f"train: {n_train:5d} ({n_train/n_total*100:.1f}%), "
          f"val: {n_val:5d} ({n_val/n_total*100:.1f}%)")

print("\n" + "=" * 80)
print(f"TRAIN: {len(train_df):6d} reactions ({len(train_df)/len(df_clean)*100:.1f}%)")
print(f"VAL:   {len(val_df):6d} reactions ({len(val_df)/len(df_clean)*100:.1f}%)")
print("=" * 80)

# Save splits
train_df[['REACTION_ID', 'REACTION_SMILES', 'EC_NUMBER']].to_csv(
    'train_split_stratified.tsv', sep='\t', index=False)
val_df[['REACTION_ID', 'REACTION_SMILES', 'EC_NUMBER']].to_csv(
    'val_split_stratified.tsv', sep='\t', index=False)

print("\n✓ Files saved: train_split_stratified.tsv, val_split_stratified.tsv")