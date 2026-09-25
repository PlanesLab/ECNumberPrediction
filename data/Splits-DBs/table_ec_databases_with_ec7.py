import pandas as pd
import os

# Same as table_ec_databases_no_ec7.py but keeps EC class 7 (Translocases)
# in both the counts and the percentage denominator.

root_path = "/scratch/jarcagniriv/ECNumberPrediction/data/Splits-DBs"
desired_order = ["MetaNetX", "KEGG", "ECREACT"]

available_databases = {
    d for d in os.listdir(root_path)
    if os.path.isdir(os.path.join(root_path, d))
}
databases = [db for db in desired_order if db in available_databases]

enzyme_classes = {
    '1': 'Oxidoreductases',
    '2': 'Transferases',
    '3': 'Hydrolases',
    '4': 'Lyases',
    '5': 'Isomerases',
    '6': 'Ligases',
    '7': 'Translocases'
}
ordered_classes = [enzyme_classes[str(i)] for i in range(1, 8)]


def class_counts_for(file_path):
    if not os.path.exists(file_path):
        return None, 0

    df = pd.read_csv(file_path, sep='\t')
    n_original = len(df)

    if 'ec' not in df.columns:
        return None, n_original

    df = df[df['ec'].notna()]
    df['ec'] = df['ec'].astype(str)
    df = df.assign(ec=df['ec'].str.split('|')).explode('ec')
    df['EC_class'] = df['ec'].str.split('.').str[0]
    df['EC_class_name'] = df['EC_class'].map(enzyme_classes)

    class_counts = df['EC_class_name'].value_counts()
    class_counts = class_counts.reindex(ordered_classes).fillna(0).astype(int)

    return class_counts, n_original


rows = []
for db_name in databases:
    db_path = os.path.join(root_path, db_name)
    for split in ('train', 'test'):
        file_path = os.path.join(db_path, f'{split}.tsv')
        class_counts, n_original = class_counts_for(file_path)
        if class_counts is None:
            continue
        n_total = int(class_counts.sum())
        for ec_class in ordered_classes:
            count = int(class_counts[ec_class])
            pct = (count / n_total * 100) if n_total else 0.0
            rows.append({
                'database': db_name,
                'split': split,
                'EC_class': ec_class,
                'count': count,
                'percentage': round(pct, 2),
                'n_total': n_total,
                'n_original': n_original,
            })

table = pd.DataFrame(rows)
output_path = os.path.join(root_path, 'ec_databases_with_ec7_table.csv')
table.to_csv(output_path, index=False)
print(f"Table saved to {output_path}")
print(table.to_string(index=False))
