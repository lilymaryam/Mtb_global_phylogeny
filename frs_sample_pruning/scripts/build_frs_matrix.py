#!/usr/bin/env python3
import sys
import pandas as pd

output_file = sys.argv[1]
input_files = sys.argv[2:]

rows = []
for f in input_files:
    with open(f) as fh:
        for line in fh:
            pair, ratio, tot, filt, result = line.strip().split('\t')
            rows.append({'sample': pair, 'ratio': int(ratio), 'value': result})

df = pd.DataFrame(rows)
matrix = df.pivot(index='sample', columns='ratio', values='value')
matrix.to_csv(output_file, sep='\t')