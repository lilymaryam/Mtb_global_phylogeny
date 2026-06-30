import argparse
import numpy as np
import matplotlib.pyplot as plt

#THRESHOLDS = [0.05, 0.1, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]
THRESHOLDS = [0.05, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
              0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]

def load_last_column(filepath):
    values = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            try:
                values.append(float(parts[-1]))
            except ValueError:
                continue
    return np.array(values)

def main():
    parser = argparse.ArgumentParser(description='Count samples exceeding threshold values and plot results.')
    parser.add_argument('file', help='Path to input data file')
    parser.add_argument('--output-tsv', default='threshold_counts.tsv', help='Output TSV filename (default: threshold_counts.tsv)')
    parser.add_argument('--output-plot', default=None, help='Save plot to this file instead of displaying (e.g. plot.png)')
    args = parser.parse_args()

    values = load_last_column(args.file)
    total = len(values)
    print(f"Loaded {total} values.\n")

    header = ['threshold', 'n_above', 'pct_above', 'n_below_or_equal', 'pct_below_or_equal']
    rows = []
    n_above_list = []

    for t in THRESHOLDS:
        n_above = int(np.sum(values > t))
        n_below = total - n_above
        pct_above = 100.0 * n_above / total
        pct_below = 100.0 * n_below / total
        n_above_list.append(n_above)
        rows.append([f'{t}', str(n_above), f'{pct_above:.2f}', str(n_below), f'{pct_below:.2f}'])

    with open(args.output_tsv, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')

    col_widths = [max(len(h), max(len(r[i]) for r in rows)) for i, h in enumerate(header)]
    fmt = '  '.join(f'{{:<{w}}}' for w in col_widths)
    print(fmt.format(*header))
    print('  '.join('-' * w for w in col_widths))
    for row in rows:
        print(fmt.format(*row))
    print(f"\nSaved to {args.output_tsv}")

    # --- plot ---
    labels = [str(t) for t in THRESHOLDS]
    cdf = [100.0 * (total - n) / total for n in n_above_list]  # % at or below threshold

    x = np.arange(len(THRESHOLDS))
    fig, ax1 = plt.subplots(figsize=(12, 5))

    ax1.bar(x, n_above_list, color='steelblue', edgecolor='white', linewidth=0.6, label='N above threshold')
    ax1.set_xlabel('Threshold')
    ax1.set_ylabel('Number of samples above threshold', color='steelblue')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=45, ha='right')

    ax2 = ax1.twinx()
    ax2.plot(x, cdf, color='tomato', linewidth=2, marker='o', markersize=4, label='CDF (% at or below)')
    ax2.set_ylabel('Cumulative % at or below threshold', color='tomato')
    ax2.tick_params(axis='y', labelcolor='tomato')
    ax2.set_ylim(0, 105)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='center right')

    plt.title('Samples above threshold with CDF')
    plt.tight_layout()

    if args.output_plot:
        plt.savefig(args.output_plot, dpi=150)
        print(f"Plot saved to {args.output_plot}")
    else:
        plt.show()

if __name__ == '__main__':
    main()