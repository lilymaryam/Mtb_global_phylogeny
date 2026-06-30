#!/usr/bin/env python3
"""
annotate_manual_prune.py

Adds a 'filter' column to a large TSV of samples with three categories
(evaluated in priority order):
  1. manual_prune  — sample ID is in the prune list (100-sample file, no header)
  2. >0.19         — frs column > 0.19
  3. <=0.19        — everything else

Usage:
    python3 annotate_manual_prune.py \
        --prune-list  ../Supplemental/SupplementalTable3.tsv \
        --large-file  /path/to/129k_samples.tsv \
        --output      /path/to/output.tsv
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prune-list", required=True)
    parser.add_argument("--large-file", required=True)
    parser.add_argument("--output",     required=True)
    parser.add_argument("--frs-cutoff", default=0.19, type=float)
    args = parser.parse_args()

    # Load prune list — single column, no header
    prune_ids = set()
    with open(args.prune_list) as fh:
        for line in fh:
            sid = line.strip()
            if sid:
                prune_ids.add(sid)
    print(f"Loaded {len(prune_ids)} IDs from prune list.")

    matched_prune = 0
    matched_high  = 0
    matched_low   = 0

    with open(args.large_file) as fh_in, open(args.output, "w") as fh_out:
        # Header line
        header = fh_in.readline().rstrip("\n")
        cols = header.split("\t")

        try:
            sample_idx = cols.index("strain")
            frs_idx    = cols.index("frs")
        except ValueError as e:
            sys.exit(f"ERROR: {e}. Columns found: {cols}")

        fh_out.write(header + "\tfilter\n")

        for line in fh_in:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            sid = parts[sample_idx]
            frs_raw = parts[frs_idx]

            if sid in prune_ids:
                label = "manual_prune"
                matched_prune += 1
            else:
                try:
                    frs_val = float(frs_raw)
                except ValueError:
                    frs_val = 0.0

                if frs_val > args.frs_cutoff:
                    label = f">{args.frs_cutoff}"
                    matched_high += 1
                else:
                    label = f"<={args.frs_cutoff}"
                    matched_low += 1

            fh_out.write(line + "\t" + label + "\n")

    total = matched_prune + matched_high + matched_low
    print(f"Processed {total:,} samples.")
    print(f"  manual_prune : {matched_prune:,}")
    print(f"  >{args.frs_cutoff}      : {matched_high:,}")
    print(f"  <={args.frs_cutoff}     : {matched_low:,}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()