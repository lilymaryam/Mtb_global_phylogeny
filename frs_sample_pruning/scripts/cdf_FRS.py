import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse

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
                continue  # skip header lines if any
    return np.array(values)

def plot(values, bins, output):
    fig, ax1 = plt.subplots(figsize=(9, 5))

    counts, edges, patches = ax1.hist(values, bins=bins, color='steelblue',
                                       edgecolor='white', linewidth=0.6, label='Count')
    ax1.set_xlabel('Value (last column)')
    ax1.set_ylabel('Count', color='steelblue')
    ax1.tick_params(axis='y', labelcolor='steelblue')

    ax2 = ax1.twinx()
    sorted_vals = np.sort(values)
    cdf = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
    ax2.plot(sorted_vals, cdf, color='tomato', linewidth=2, label='CDF')
    ax2.set_ylabel('Cumulative proportion', color='tomato')
    ax2.tick_params(axis='y', labelcolor='tomato')
    ax2.set_ylim(0, 1.05)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

    plt.title('Histogram and CDF of last column')
    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=150)
        print(f"Saved to {output}")
    else:
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot histogram + CDF of the last column in a whitespace-delimited file.')
    parser.add_argument('--file', help='Path to input data file')
    parser.add_argument('--bins', type=int, default=20, help='Number of histogram bins (default: 20)')
    parser.add_argument('--output', default=None, help='Save plot to this file instead of displaying (e.g. plot.png)')
    args = parser.parse_args()

    values = load_last_column(args.file)
    print(f"Loaded {len(values)} values. Min={values.min():.4f}, Max={values.max():.4f}, Mean={values.mean():.4f}")
    plot(values, bins=args.bins, output=args.output)