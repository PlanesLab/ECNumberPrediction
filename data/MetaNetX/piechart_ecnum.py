import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import colorsys
import matplotlib.colors as mcolors  

csv_path = '/scratch/jarcagniriv/Case2/MetaNetX/train_reactions.tsv'

if os.path.exists(csv_path):
    df = pd.read_csv(csv_path, sep='\t')

    if 'EC_number' in df.columns:
        df = df[df['EC_number'].notna()]
        
        # Split EC numbers and expand DataFrame
        df['EC_number'] = df['EC_number'].astype(str)
        df = df.assign(EC_number=df['EC_number'].str.split('|')).explode('EC_number')
        df['EC_class'] = df['EC_number'].str.split('.').str[0]
        df['EC_subclass'] = df['EC_number'].apply(lambda x: '.'.join(x.split('.')[:2]) if '.' in x else None)

        enzyme_classes = {
            '1': 'Oxidoreductases',
            '2': 'Transferases',
            '3': 'Hydrolases',
            '4': 'Lyases',
            '5': 'Isomerases',
            '6': 'Ligases',
            '7': 'Translocases'
        }
        df['EC_class_name'] = df['EC_class'].map(enzyme_classes)

        # Count occurrences
        class_counts = df['EC_class_name'].value_counts()
        ordered_classes = [enzyme_classes[str(i)] for i in range(1, 8)]
        class_counts = class_counts.reindex(ordered_classes).dropna()
        subclass_counts = df.groupby(['EC_class_name', 'EC_subclass']).size()

        class_color_map = {
            'Oxidoreductases': '#f3aa7c',
            'Transferases':    '#a5dbbc',
            'Hydrolases':      '#8adadf',
            'Lyases':          '#f3e57c',
            'Isomerases':      '#E58A98',
            'Ligases':         '#a5a9db',
            'Translocases':    '#C2C2C2'
        }

        def adjust_color_brightness(color, factor):
            r, g, b = color[:3]
            h, l, s = colorsys.rgb_to_hls(r, g, b)
            new_l = min(max(l * factor, 0.2), 0.95) 
            return colorsys.hls_to_rgb(h, new_l, s)

        inner_labels = class_counts.index.tolist()
        inner_sizes = class_counts.values
        inner_colors = [class_color_map[class_name] for class_name in inner_labels]

        outer_labels = []
        outer_sizes = []
        outer_colors = []

        for class_name in inner_labels:
            if class_name in subclass_counts.index.levels[0]:
                class_subclasses = subclass_counts[class_name]
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
        threshold_pct = 1
        outer_labels_adjusted = [
            label if (size / total_outer * 100) >= threshold_pct else ''
            for label, size in zip(outer_labels, outer_sizes)
        ]

        def inner_autopct(pct):
            return ('%1.1f%%' % pct) if pct >= 1 else ''

        fig, ax = plt.subplots(figsize=(10, 10))

        ax.pie(
            inner_sizes,
            autopct=inner_autopct,
            startangle=90,
            colors=inner_colors,
            wedgeprops=dict(width=0.4, edgecolor='white'),
            radius=0.7,
            pctdistance=0.85,
            textprops={'fontsize': 10}
        )

        ax.pie(
            outer_sizes,
            labels=outer_labels_adjusted,
            startangle=90,
            colors=outer_colors,
            radius=1.0,
            wedgeprops=dict(width=0.25, edgecolor='white'),
            textprops={'fontsize': 8},
            labeldistance=1.02
        )

        ax.legend(inner_labels, title="Enzyme Class", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
        plt.title("MetaNetX Test EC Number Distribution")
        plt.tight_layout()
        plt.savefig('/scratch/jarcagniriv/Case2/MetaNetX/ecdist2025_train.png', dpi=600, bbox_inches='tight')

    else:
        print("The 'EC Number' column does not exist in the dataset.")
else:
    print(f"The file {csv_path} does not exist.")
