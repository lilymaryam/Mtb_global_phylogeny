#!/usr/bin/env python3
import sys
import json

r1, r2, sample_name, output_file = sys.argv[1:5]

inputs = {
    "myco.paired_fastq_sets": [[r1, r2]],
    "myco.output_sample_name": sample_name
}

with open(output_file, 'w') as f:
    json.dump(inputs, f, indent=2)