import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import colorsys
import matplotlib.colors as mcolors  

# ---------- CONFIG ----------
train_csv = '/scratch/jarcagniriv/Case2/MetaNetX/train_reactions_expanded.csv'
test_csv  = '/scratch/jarcagniriv/Case2/MetaNetX/test_reactions.csv'
output_png = '/scratch/jarcagniriv/Case2/MetaNetX/ecdist2025_train_test.jpg'

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

# ---------- FUNCTIONS ----------
def adjust_color_brightness(color, factor):
    r, g, b = color[:3]
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    new_l = min(max(l * factor, 0.2), 0.95) 
    return colorsys.hls_to_rgb(h, new_l, s)

def process_file(csv_path):
    if not os.path.exists(csv_path):
        print(f"File {csv_path} not found.")
        return None
    
    df = pd.read_csv(csv_path)
    if 'EC_number' not in df.columns:
        print(f"No EC_number column in {csv_path}.")
        return None

    df = df[df['EC_number'].notna()]
    df['EC_number'] = df['EC_number'].astype(str)
    df = df.assign(EC_number=df['EC_number'].str.split('|')).explode('EC_number')
    df['EC_class'] = df['EC_number'].str.split('.').str[0]
    df['EC_subclass'] = df['EC_number'].apply(lambda x: '.'.join(x.split('.')[:2]) if '.' in x else None)
    df['EC_class_name'] = df['EC_class'].map(enzyme_classes)

    class_counts = df['EC_class_name'].value_counts()
    ordered_classes = [enzyme_classes[str(i)] for i in range(1, 8)]
    class_counts = class_counts.reindex(ordered_classes).dropna()
    subclass_counts = df.groupby(['EC_class_name', 'EC_subclass']).size()

    return class_counts, subclass_counts

def plot_pie(ax, class_counts, subclass_counts, title, label_letter):
    inner_labels = class_counts.index.tolist()
    inner_sizes = class_counts.values
    inner_colors = [class_color_map[class_name] for class_name in inner_labels]

    outer_labels, outer_sizes, outer_colors = [], [], []
    for class_name in inner_labels:
        if class_name in subclass_counts.index.levels[0]:
            class_subclasses = subclass_counts[class_name]
            max_count = class_subclasses.max()
            min_count = class_subclasses.min()
            range_count = max_count - min_count if max_count != min_count else 1

            base_color_rgb = mcolors.to_rgb(class_color_map[class_name])
            for subclass, count in class_subclasses.items():
                intensity = 0.95 - ((count - min_count) / range_count) * 0.3
                shade = adjust_color_brightness(base_color_rgb, intensity)
                outer_labels.append(subclass)
                outer_sizes.append(count)
                outer_colors.append(shade)

    total_outer = sum(outer_sizes)
    threshold_pct = 2
    outer_labels_adjusted = [
        label if (size / total_outer * 100) >= threshold_pct else ''
        for label, size in zip(outer_labels, outer_sizes)
    ]

    # Track wedges + hidden percentages
    hidden_labels = []

    def inner_autopct(pct):
        if pct >= 2.9:
            return f'{pct:.0f}%'
        else:
            hidden_labels.append(pct)
            return ''  # don’t show inside pie

    wedges, texts, autotexts = ax.pie(
        inner_sizes,
        autopct=inner_autopct,
        startangle=90,
        colors=inner_colors,
        wedgeprops=dict(width=0.35, edgecolor='white'),
        radius=0.7,
        pctdistance=0.8,  
        textprops={'fontsize': 20, 'weight': 'normal'}
    )

    ax.pie(
        outer_sizes,
        labels=outer_labels_adjusted,
        startangle=90,
        colors=outer_colors,
        radius=1.0,
        wedgeprops=dict(width=0.25, edgecolor='white'),
        textprops={'fontsize': 15},
        labeldistance=1.02
    )

    # 🔥 Add arrows for hidden small percentages
    for i, w in enumerate(wedges):
        pct = inner_sizes[i] / sum(inner_sizes) * 100
        if pct < 2.9:
            ang = (w.theta2 + w.theta1) / 2.0  # angle of wedge center
            x = np.cos(np.deg2rad(ang)) * 0.85
            y = np.sin(np.deg2rad(ang)) * 0.85
            start_factor = 1.2  # start further from the wedge
            line_factor = 1.31  # end of the line
            text_factor = 1.36   # text position

            ax.plot([start_factor*x, line_factor*x], [start_factor*y, line_factor*y], lw=1.2, color='black')
            ax.text(text_factor*x, text_factor*y, f"{pct:.1f}%", ha='center', va='center', fontsize=20)

    ax.set_title(title, fontsize=26, weight='normal')
    ax.text(-0.2, 1.15, label_letter, fontsize=30, weight='bold', transform=ax.transAxes)

# ---------- MAIN ----------
train_data = process_file(train_csv)
test_data = process_file(test_csv)

if train_data and test_data:
    fig, axes = plt.subplots(1, 2, figsize=(22, 12))

    plot_pie(axes[0], train_data[0], train_data[1], "Train EC Number Distribution", "A")
    plot_pie(axes[1], test_data[0], test_data[1], "Test EC Number Distribution", "B")

    # Legend with no title, very big
    axes[1].legend(
        train_data[0].index,
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1),
        fontsize=25 
    )

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.show()
