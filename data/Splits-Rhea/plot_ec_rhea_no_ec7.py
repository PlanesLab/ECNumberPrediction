import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import colorsys
import matplotlib.colors as mcolors

# ================================
# ROOT DATABASE DIRECTORY
# ================================
root_path = "/scratch/jarcagniriv/ECNumberPrediction/data/Splits-Rhea"

# Automatically detect split folders
# Desired order
desired_order = ["Stratified", "Time", "Scaffold"]

# Detect existing folders
available_databases = {
    d for d in os.listdir(root_path)
    if os.path.isdir(os.path.join(root_path, d))
}

# Keep only those in desired order AND present
databases = [db for db in desired_order if db in available_databases]

print(f"Using databases in order: {databases}")


print(f"Detected databases: {databases}")

# ================================
# Enzyme class mapping
# ================================
enzyme_classes = {
    '1': 'Oxidoreductases',
    '2': 'Transferases',
    '3': 'Hydrolases',
    '4': 'Lyases',
    '5': 'Isomerases',
    '6': 'Ligases',
    '7': 'Translocases'
}

class_color_map = {
    'Oxidoreductases': '#f3aa7c',
    'Transferases':    '#a5dbbc',
    'Hydrolases':      '#8adadf',
    'Lyases':          '#f3e57c',
    'Isomerases':      '#E58A98',
    'Ligases':         '#a5a9db',
    'Translocases':    '#C2C2C2'
}

# ================================
# Utility Functions
# ================================

def adjust_color_brightness(color, factor):
    r, g, b = color[:3]
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    new_l = min(max(l * factor, 0.2), 0.95)
    return colorsys.hls_to_rgb(h, new_l, s)


def class_and_subclass_counts(df):
    """From a raw train/test dataframe, return (class_counts, subclass_counts)
    Series (class 7 / Translocases excluded), keyed the same way regardless
    of which seed/file the rows came from."""
    ec_col = 'ec' if 'ec' in df.columns else ('EC_NUMBER' if 'EC_NUMBER' in df.columns else None)
    if ec_col is None:
        return None, None
    df = df.rename(columns={ec_col: 'ec'})

    df = df[df['ec'].notna()]
    df['ec'] = df['ec'].astype(str)
    df = df.assign(ec=df['ec'].str.split('|')).explode('ec')
    df['EC_class'] = df['ec'].str.split('.').str[0]
    df = df[df['EC_class'] != '7']  # exclude class 7 (Translocases); percentages recompute over the remaining classes
    df['EC_subclass'] = df['ec'].apply(lambda x: '.'.join(x.split('.')[:2]) if '.' in x else None)
    df['EC_class_name'] = df['EC_class'].map(enzyme_classes)

    class_counts = df['EC_class_name'].value_counts()
    ordered_classes = [enzyme_classes[str(i)] for i in range(1, 7)]
    class_counts = class_counts.reindex(ordered_classes).fillna(0)
    subclass_counts = df.groupby(['EC_class_name', 'EC_subclass']).size()

    return class_counts, subclass_counts


def build_wedge_data(class_counts, subclass_counts, total_n):
    """Turn (class_counts, subclass_counts) -- either raw counts from one
    file, or seed-averaged fractional 'counts' -- into the inner/outer wedge
    arrays the plotting code expects. Scale-invariant: works the same whether
    the values are integer counts or averaged fractions."""
    class_counts = class_counts[class_counts > 0]

    inner_labels = class_counts.index.tolist()
    inner_sizes = class_counts.values
    inner_colors = [class_color_map[class_name] for class_name in inner_labels]

    outer_labels = []
    outer_sizes = []
    outer_colors = []

    for class_name in inner_labels:
        if class_name in subclass_counts.index.levels[0]:
            class_subclasses = subclass_counts[class_name]
            class_subclasses = class_subclasses[class_subclasses > 0]
            max_count = class_subclasses.max()
            min_count = class_subclasses.min()
            range_count = max_count - min_count if max_count != min_count else 1

            base_color_hex = class_color_map[class_name]
            base_color_rgb = mcolors.to_rgb(base_color_hex)

            for subclass, count in class_subclasses.items():
                intensity = 0.95 - ((count - min_count) / range_count) * 0.3
                shade = adjust_color_brightness(base_color_rgb, intensity)
                outer_labels.append(subclass)
                outer_sizes.append(count)
                outer_colors.append(shade)

    total_outer = sum(outer_sizes)
    threshold_pct = 1.5
    outer_labels_adjusted = [
        label if (size / total_outer * 100) >= threshold_pct else ''
        for label, size in zip(outer_labels, outer_sizes)
    ]

    return inner_labels, inner_sizes, inner_colors, outer_labels_adjusted, outer_sizes, outer_colors, total_n


def process_ec_data(file_path):
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None, None, None, None, None, None, 0

    df = pd.read_csv(file_path, sep='\t')
    total_n = len(df)

    class_counts, subclass_counts = class_and_subclass_counts(df)
    if class_counts is None:
        print(f"ec column not found in {file_path}")
        return None, None, None, None, None, None, total_n

    return build_wedge_data(class_counts, subclass_counts, total_n)


def process_ec_data_avg_seeds(seed_dirs, split_name):
    """Average EC class/subclass composition across several seed dirs (each
    holding its own train.tsv/test.tsv), for split strategies (Stratified,
    Scaffold) that were re-run with multiple random seeds. Each seed
    contributes its class/subclass proportions (not raw counts) with equal
    weight, then those proportions are averaged -- so a seed with a
    slightly larger/smaller split doesn't get over/under-weighted."""
    class_frac_sum = None
    subclass_frac_sum = None
    totals = []

    for seed_dir in seed_dirs:
        file_path = os.path.join(seed_dir, f'{split_name}.tsv')
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue
        df = pd.read_csv(file_path, sep='\t')
        raw_total = len(df)  # original database row count (incl. class 7), matching process_ec_data's n
        class_counts, subclass_counts = class_and_subclass_counts(df)
        if class_counts is None:
            print(f"ec column not found in {file_path}")
            continue

        totals.append(raw_total)
        class_frac = class_counts / class_counts.sum()
        subclass_frac = subclass_counts / class_counts.sum()

        class_frac_sum = class_frac if class_frac_sum is None else class_frac_sum.add(class_frac, fill_value=0)
        subclass_frac_sum = subclass_frac if subclass_frac_sum is None else subclass_frac_sum.add(subclass_frac, fill_value=0)

    if not totals:
        return None, None, None, None, None, None, 0

    n_seeds = len(totals)
    class_frac_avg = class_frac_sum / n_seeds
    subclass_frac_avg = subclass_frac_sum / n_seeds
    avg_n = round(sum(totals) / n_seeds)

    return build_wedge_data(class_frac_avg, subclass_frac_avg, avg_n)


def add_external_labels_for_small_wedges(ax, wedges, sizes, threshold=6.0, cluster_gap=25):
    total = sum(sizes)

    # Collect small wedges
    small = []
    for wedge, size in zip(wedges, sizes):
        pct = 100 * size / total
        if pct < threshold:
            ang = (wedge.theta2 - wedge.theta1) / 2. + wedge.theta1
            small.append((wedge, pct, ang))

    # Sort by angle so neighbors get different radii
    small.sort(key=lambda x: x[2])

    # Group angularly-close labels into clusters so isolated small wedges
    # stay near the pie, while crowded neighbors (e.g. adjacent Isomerases/
    # Ligases) get pushed apart radially AND nudged apart in angle so their
    # text doesn't collide.
    clusters = []
    for wedge, pct, ang in small:
        if clusters and (ang - clusters[-1][-1][2]) < cluster_gap:
            clusters[-1].append((wedge, pct, ang))
        else:
            clusters.append([(wedge, pct, ang)])

    for cluster in clusters:
        n = len(cluster)
        for i, (wedge, pct, ang) in enumerate(cluster):
            # Spread members of a crowded cluster apart in angle too, so
            # their radial lines (and text) don't sit on top of each other.
            if n > 1:
                spread = 7  # degrees between adjacent cluster members
                ang_offset = (i - (n - 1) / 2) * spread
            else:
                ang_offset = 0
            label_ang = ang + ang_offset

            # Keep exact angular position for the connector's start point
            x0 = np.cos(np.deg2rad(ang))
            y0 = np.sin(np.deg2rad(ang))
            inner_radius = wedge.r - wedge.width / 2
            x_start = x0 * inner_radius
            y_start = y0 * inner_radius

            # Label/connector end point uses the nudged angle plus a radius
            # that grows with position in the cluster.
            x = np.cos(np.deg2rad(label_ang))
            y = np.sin(np.deg2rad(label_ang))
            label_distance = 1.10 + (i * 0.08)
            x_end = x * label_distance
            y_end = y * label_distance

            ha = 'left' if x >= 0 else 'right'

            # Straight connector
            ax.plot([x_start, x_end],
                    [y_start, y_end],
                    color='gray',
                    lw=0.8)

            ax.text(x_end, y_end,
                    f'{pct:.1f}%',
                    ha=ha,
                    va='center',
                    fontsize=12,
                    weight='bold')



def inner_autopct(pct):
    return ('%1.1f%%' % pct) if pct >= 6.0 else ''


def add_outer_subclass_labels(ax, wedges, labels, label_distance=1.04, spread=5.5, cluster_gap=8):
    """Place the outer ring's subclass code labels (e.g. '5.1', '5.3') right
    next to the ring, like matplotlib's built-in labels=/labeldistance=, but
    nudge angularly-crowded neighbors apart in angle (no radial push-out, no
    leader lines) so adjacent thin subclass wedges don't overlap."""
    items = []
    for wedge, label in zip(wedges, labels):
        if not label:
            continue
        ang = (wedge.theta2 - wedge.theta1) / 2. + wedge.theta1
        items.append((wedge, label, ang))

    items.sort(key=lambda x: x[2])

    clusters = []
    for wedge, label, ang in items:
        if clusters and (ang - clusters[-1][-1][2]) < cluster_gap:
            clusters[-1].append((wedge, label, ang))
        else:
            clusters.append([(wedge, label, ang)])

    for cluster in clusters:
        n = len(cluster)
        for i, (wedge, label, ang) in enumerate(cluster):
            ang_offset = (i - (n - 1) / 2) * spread if n > 1 else 0
            label_ang = ang + ang_offset

            x = np.cos(np.deg2rad(label_ang))
            y = np.sin(np.deg2rad(label_ang))
            x_end = x * label_distance
            y_end = y * label_distance

            ha = 'left' if x >= 0 else 'right'

            ax.text(x_end, y_end, label, ha=ha, va='center', fontsize=11)


# ================================
# CREATE FIGURE
# ================================

n_databases = len(databases)
fig, axes = plt.subplots(n_databases, 2, figsize=(13, 5 * n_databases))

if n_databases == 1:
    axes = np.array([axes])  # ensure 2D structure

for row_idx, db_name in enumerate(databases):

    db_path = os.path.join(root_path, db_name)

    train_path = os.path.join(db_path, 'train.tsv')
    test_path = os.path.join(db_path, 'test.tsv')

    # Splits re-run with multiple random seeds (Stratified, Scaffold) get
    # averaged across all their seed_splits/seed* dirs; Time is a single
    # deterministic chronological split (only ever has seed0), so it just
    # uses its one file as before.
    seed_splits_dir = os.path.join(db_path, 'seed_splits')
    seed_dirs = sorted(
        os.path.join(seed_splits_dir, d) for d in os.listdir(seed_splits_dir)
    ) if os.path.isdir(seed_splits_dir) else []
    is_multi_seed = len(seed_dirs) > 1
    n_seeds = len(seed_dirs)

    # ---- TRAIN ----
    if is_multi_seed:
        train_data = process_ec_data_avg_seeds(seed_dirs, 'train')
    else:
        train_data = process_ec_data(train_path)
    ax_train = axes[row_idx, 0]

    if train_data[0] is not None:
        inner_labels, inner_sizes, inner_colors, outer_labels, outer_sizes, outer_colors, n_train = train_data

        wedges, _, _ = ax_train.pie(
            inner_sizes,
            autopct=inner_autopct,
            startangle=90,
            colors=inner_colors,
            wedgeprops=dict(width=0.4, edgecolor='white', linewidth=1.5),
            radius=0.7,
            pctdistance=0.88,
            textprops={'fontsize': 12, 'weight': 'bold'}
        )

        add_external_labels_for_small_wedges(ax_train, wedges, inner_sizes)

        outer_wedges, _ = ax_train.pie(
            outer_sizes,
            startangle=90,
            colors=outer_colors,
            radius=1.0,
            wedgeprops=dict(width=0.25, edgecolor='white', linewidth=1)
        )
        add_outer_subclass_labels(ax_train, outer_wedges, outer_labels)

        train_title_n = f'n≈{n_train:,}' if is_multi_seed else f'n={n_train:,}'
        ax_train.set_title(f'{db_name} - Train ({train_title_n})', fontsize=17, fontweight='bold', pad=28)

    # ---- TEST ----
    if is_multi_seed:
        test_data = process_ec_data_avg_seeds(seed_dirs, 'test')
    else:
        test_data = process_ec_data(test_path)
    ax_test = axes[row_idx, 1]

    if test_data[0] is not None:
        inner_labels, inner_sizes, inner_colors, outer_labels, outer_sizes, outer_colors, n_test = test_data

        wedges, _, _ = ax_test.pie(
            inner_sizes,
            autopct=inner_autopct,
            startangle=90,
            colors=inner_colors,
            wedgeprops=dict(width=0.4, edgecolor='white', linewidth=1.5),
            radius=0.7,
            pctdistance=0.88,
            textprops={'fontsize': 12, 'weight': 'bold'}
        )

        add_external_labels_for_small_wedges(ax_test, wedges, inner_sizes)

        outer_wedges, _ = ax_test.pie(
            outer_sizes,
            startangle=90,
            colors=outer_colors,
            radius=1.0,
            wedgeprops=dict(width=0.25, edgecolor='white', linewidth=1)
        )
        add_outer_subclass_labels(ax_test, outer_wedges, outer_labels)

        test_title_n = f'n≈{n_test:,}' if is_multi_seed else f'n={n_test:,}'
        ax_test.set_title(f'{db_name} - Test ({test_title_n})', fontsize=17, fontweight='bold', pad=28)


# ================================
# LEGEND
# ================================

ordered_enzyme_classes = [enzyme_classes[str(i)] for i in range(1, 7)]
legend_colors = [class_color_map[cls] for cls in ordered_enzyme_classes]
legend_handles = [plt.Rectangle((0, 0), 1, 1, fc=color) for color in legend_colors]

legend = fig.legend(
    legend_handles,
    ordered_enzyme_classes,
    loc='center right',
    bbox_to_anchor=(0.98, 0.5),
    fontsize=17,
    title_fontsize=18
)

plt.tight_layout(rect=[0, 0, 0.77, 1.0])
output_path = os.path.join(root_path, "ec_rhea_no_ec7.png")
plt.savefig(output_path, dpi=600, bbox_inches='tight', pad_inches=0.6)

print(f"Figure saved to {output_path}")
