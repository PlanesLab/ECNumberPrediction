import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import numpy as np

def get_bemis_murcko_scaffold(smiles):
    """Extract Bemis-Murcko scaffold from a SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except:
        return None

def extract_product_smiles(reaction_smiles):
    """Extract product from reaction SMILES (after >>)"""
    try:
        parts = reaction_smiles.split('>>')
        if len(parts) == 2:
            products = parts[1].split('.')
            return products[0].strip()
        return None
    except:
        return None

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

# Extract products and scaffolds
print("\nExtracting Bemis-Murcko scaffolds...")
df['product_smiles'] = df['REACTION_SMILES'].apply(extract_product_smiles)
df['scaffold'] = df['product_smiles'].apply(get_bemis_murcko_scaffold)
df['ec_subsubclass'] = df['EC_NUMBER'].apply(get_ec_subsubclass)

# Remove rows where scaffold extraction failed
df_clean = df.dropna(subset=['scaffold', 'ec_subsubclass'])
print(f"Reactions with valid scaffolds: {len(df_clean)}")

# Group by EC subsubclass and perform split
print("\nPerforming scaffold-aware split (90% train, 10% val)...")
print("-" * 80)

train_indices = []
val_indices = []

ec_groups = df_clean.groupby('ec_subsubclass')

for ec_class, group in ec_groups:
    # Get unique scaffolds in this EC subsubclass
    scaffolds = group['scaffold'].unique()
    n_scaffolds = len(scaffolds)
    
    # Calculate 90-10 split
    n_train = max(1, int(0.9 * n_scaffolds))
    
    # Shuffle scaffolds
    np.random.seed(42)
    shuffled_scaffolds = np.random.permutation(scaffolds)
    
    train_scaffolds = set(shuffled_scaffolds[:n_train])
    val_scaffolds = set(shuffled_scaffolds[n_train:])
    
    # Assign reactions based on scaffold
    group_train = group[group['scaffold'].isin(train_scaffolds)]
    group_val = group[group['scaffold'].isin(val_scaffolds)]
    
    train_indices.extend(group_train.index.tolist())
    val_indices.extend(group_val.index.tolist())
    
    print(f"EC {ec_class}: {n_scaffolds:4d} scaffolds → "
          f"train: {len(train_scaffolds):4d} ({len(group_train):5d} rxns), "
          f"val: {len(val_scaffolds):4d} ({len(group_val):5d} rxns)")

# Create train and val datasets
train_df = df_clean.loc[train_indices]
val_df = df_clean.loc[val_indices]

print("\n" + "=" * 80)
print(f"TRAIN: {len(train_df):6d} reactions ({len(train_df)/len(df_clean)*100:.1f}%)")
print(f"VAL:   {len(val_df):6d} reactions ({len(val_df)/len(df_clean)*100:.1f}%)")
print("=" * 80)

# Save splits
train_df[['REACTION_ID', 'REACTION_SMILES', 'EC_NUMBER']].to_csv(
    'train_split.tsv', sep='\t', index=False)
val_df[['REACTION_ID', 'REACTION_SMILES', 'EC_NUMBER']].to_csv(
    'val_split.tsv', sep='\t', index=False)

# Statistics
train_scaffolds = set(train_df['scaffold'])
val_scaffolds = set(val_df['scaffold'])
overlap = train_scaffolds & val_scaffolds

print(f"\nUnique scaffolds in train: {len(train_scaffolds)}")
print(f"Unique scaffolds in val:   {len(val_scaffolds)}")
print(f"Scaffold overlap:          {len(overlap)}")
print("\n✓ Files saved: train_split.tsv, val_split.tsv")