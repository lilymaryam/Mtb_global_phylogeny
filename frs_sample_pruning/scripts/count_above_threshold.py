#!/usr/bin/env python3

import argparse


def main():
    parser = argparse.ArgumentParser(description='Count samples with FRS above each ratio threshold')
    parser.add_argument('-s', '--samples', required=True, help='Sample file (sample_id, frs, lineage)')
    parser.add_argument('-t', '--thresholds', required=True, help='Summary TSV from FRS analysis (ratio, mean, stddev, n)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    args = parser.parse_args()
    
    # Load sample FRS values
    sample_frs = []
    with open(args.samples) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 4:
                sample_frs.append(float(parts[3]))
    
    print(f"Loaded {len(sample_frs)} sample FRS values")
    
    # Load thresholds (mean FRS per ratio)
    thresholds = []
    with open(args.thresholds) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                ratio = parts[0]
                mean_frs = float(parts[1])
                thresholds.append((ratio, mean_frs))
    
    with open(args.output, 'w') as out:
        out.write("ratio\tthreshold\tn_above\ttotal\tpct_above\n")
        for ratio, threshold in thresholds:
            n_above = sum(1 for frs in sample_frs if frs > threshold)
            pct = n_above / len(sample_frs) * 100
            out.write(f"{ratio}\t{threshold:.4f}\t{n_above}\t{len(sample_frs)}\t{pct:.2f}\n")
    
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()