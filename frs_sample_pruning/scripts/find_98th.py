#!/usr/bin/env python3

import argparse
import os
import pandas as pd


def existing_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path


def main():
    parser = argparse.ArgumentParser(
        description='Calculate 98th percentile threshold per ratio column'
    )
    parser.add_argument('-i', '--input', required=True, type=existing_file,
                        help='FRS matrix TSV file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output thresholds TSV file')
    parser.add_argument('-p', '--percentile', type=float, default=0.98,
                        help='Percentile threshold (default: 0.98)')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t', index_col=0, na_values='NA')
    df = df.apply(pd.to_numeric, errors='coerce')

    if df.empty:
        raise ValueError("Input matrix is empty or has no numeric data")

    stats = pd.DataFrame({
        'p98': df.quantile(args.percentile),
        'n': df.count(),
    })

    stats.to_csv(args.output, sep='\t')
    print(f"Wrote thresholds to {args.output}")


if __name__ == "__main__":
    main()