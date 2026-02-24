import pandas as pd

# Load old database (train)
print("Loading old database (2018)...")
df_old = pd.read_csv('/scratch/jarcagniriv/Case1-Rhea/Database/oldDB-2018-11-26/smiles_vs_ec_2018_DIRECTION_AWARE.tsv', sep='\t')
print(f"Old DB reactions: {len(df_old)}")

# Load new database (full)
print("\nLoading new database...")
df_new = pd.read_csv('/scratch/jarcagniriv/Case1-Rhea/Database/smiles_vs_ec_with_id.tsv', sep='\t')
print(f"New DB reactions: {len(df_new)}")

# Get REACTION_IDs from old database
old_ids = set(df_old['REACTION_ID'])
print(f"\nUnique REACTION_IDs in old DB: {len(old_ids)}")

# Split new database into train (old IDs) and test (new IDs)
train_df = df_new[df_new['REACTION_ID'].isin(old_ids)]
test_df = df_new[~df_new['REACTION_ID'].isin(old_ids)]

print("\n" + "=" * 80)
print("TEMPORAL SPLIT RESULTS:")
print("-" * 80)
print(f"TRAIN (old reactions):  {len(train_df):6d} reactions ({len(train_df)/len(df_new)*100:.1f}%)")
print(f"TEST  (new reactions):  {len(test_df):6d} reactions ({len(test_df)/len(df_new)*100:.1f}%)")
print("=" * 80)

# Save splits
train_df.to_csv('train_temporal.tsv', sep='\t', index=False)
test_df.to_csv('test_temporal.tsv', sep='\t', index=False)

print("\n✓ Files saved:")
print("  - train_temporal.tsv (reactions from 2018 DB)")
print("  - test_temporal.tsv (new reactions added after 2018)")

# Additional statistics
print("\n" + "=" * 80)
print("ADDITIONAL INFO:")
print(f"Reactions in old DB not found in new DB: {len(old_ids - set(train_df['REACTION_ID']))}")
print("=" * 80)