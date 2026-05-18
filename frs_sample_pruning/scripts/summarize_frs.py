#!/usr/bin/env python3

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Calculate mean and stddev per ratio column')
    parser.add_argument('matrix', help='FRS matrix TSV file')
    parser.add_argument('output', help='Output summary TSV file')
    args = parser.parse_args()
    
    df = pd.read_csv(args.matrix, sep='\t', index_col=0, na_values='NA')
    df = df.apply(pd.to_numeric, errors='coerce')
    
    stats = pd.DataFrame({
        'mean': df.mean(),
        'stddev': df.std(),
        'n': df.count(),
    })
    
    stats.to_csv(args.output, sep='\t')
    print(f"Wrote summary to {args.output}")


if __name__ == "__main__":
    main()