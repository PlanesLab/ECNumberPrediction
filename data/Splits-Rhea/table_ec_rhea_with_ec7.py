import pandas as pd
import os

# Same as table_ec_rhea_no_ec7.py but keeps EC class 7 (Translocases) in
# both the counts and the percentage denominator. Stratified and Scaffold
# are seed-averaged across their 3 seed_splits/seed* dirs; Time uses its
# single file.

root_path = "/scratch/jarcagniriv/ECNumberPrediction/data/Splits-Rhea"
desired_order = ["Stratified", "Time", "Scaffold"]

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


def class_counts_from_df(df):
    ec_col = 'ec' if 'ec' in df.columns else ('EC_NUMBER' if 'EC_NUMBER' in df.columns else None)
    if ec_col is None:
        return None
    df = df.rename(columns={ec_col: 'ec'})
    df = df[df['ec'].notna()]
    df['ec'] = df['ec'].astype(str)
    df = df.assign(ec=df['ec'].str.split('|')).explode('ec')
    df['EC_class'] = df['ec'].str.split('.').str[0]
    df['EC_class_name'] = df['EC_class'].map(enzyme_classes)

    class_counts = df['EC_class_name'].value_counts()
    return class_counts.reindex(ordered_classes).fillna(0)


def single_file(file_path):
    if not os.path.exists(file_path):
        return None, 0
    df = pd.read_csv(file_path, sep='\t')
    n_original = len(df)
    class_counts = class_counts_from_df(df)
    if class_counts is None:
        return None, n_original
    return class_counts.astype(int), n_original


def seed_averaged(seed_dirs, split_name):
    class_frac_sum = None
    totals = []
    for seed_dir in seed_dirs:
        file_path = os.path.join(seed_dir, f'{split_name}.tsv')
        if not os.path.exists(file_path):
            continue
        df = pd.read_csv(file_path, sep='\t')
        totals.append(len(df))
        class_counts = class_counts_from_df(df)
        if class_counts is None:
            continue
        class_frac = class_counts / class_counts.sum()
        class_frac_sum = class_frac if class_frac_sum is None else class_frac_sum.add(class_frac, fill_value=0)

    if not totals:
        return None, 0, 0
    n_seeds = len(totals)
    class_frac_avg = class_frac_sum / n_seeds
    avg_n_original = round(sum(totals) / n_seeds)
    return class_frac_avg, avg_n_original, n_seeds


rows = []
for db_name in databases:
    db_path = os.path.join(root_path, db_name)
    train_path = os.path.join(db_path, 'train.tsv')
    test_path = os.path.join(db_path, 'test.tsv')

    seed_splits_dir = os.path.join(db_path, 'seed_splits')
    seed_dirs = sorted(
        os.path.join(seed_splits_dir, d) for d in os.listdir(seed_splits_dir)
    ) if os.path.isdir(seed_splits_dir) else []
    is_multi_seed = len(seed_dirs) > 1
    n_seeds = len(seed_dirs) if is_multi_seed else 1

    for split, file_path in (('train', train_path), ('test', test_path)):
        if is_multi_seed:
            class_data, avg_n_original, n_seeds = seed_averaged(seed_dirs, split)
            is_avg = True
        else:
            class_data, avg_n_original = single_file(file_path)
            is_avg = False

        if class_data is None:
            continue

        n_total = float(class_data.sum()) if is_avg else int(class_data.sum())

        for ec_class in ordered_classes:
            value = class_data[ec_class]
            if is_avg:
                pct = value / n_total * 100 if n_total else 0.0
                count_display = round(value * avg_n_original / n_total) if n_total else 0
            else:
                pct = (value / n_total * 100) if n_total else 0.0
                count_display = int(value)

            rows.append({
                'split_strategy': db_name,
                'dataset': split,
                'EC_class': ec_class,
                'count': count_display,
                'percentage': round(pct, 2),
                'n_original': avg_n_original,
                'n_seeds': n_seeds,
                'seed_averaged': is_avg,
            })

table = pd.DataFrame(rows)
output_path = os.path.join(root_path, 'ec_rhea_with_ec7_table.csv')
table.to_csv(output_path, index=False)
print(f"Table saved to {output_path}")
print(table.to_string(index=False))
