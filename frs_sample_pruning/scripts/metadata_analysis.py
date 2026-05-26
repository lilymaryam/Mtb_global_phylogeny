#!/usr/bin/env python3
"""
Threshold analysis of FRS metadata file.

For each threshold column, reports:
  - Total samples cut (FRS > threshold)
  - Mean terminal branch length of cut samples
  - Median terminal branch length of cut samples
  - % of cut samples with terminal_branch_length <= 2
"""

import sys
import csv
import statistics

def parse_metadata(filepath):
    """Parse the tab-separated metadata file."""
    rows = []
    thresholds = []

    with open(filepath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames

        # Collect threshold column names (everything between 'frs' and 'terminal_branch_length')
        frs_idx = fieldnames.index("frs")
        tbl_idx = fieldnames.index("terminal_branch_length")
        thresholds = fieldnames[frs_idx + 1 : tbl_idx]

        for row in reader:
            rows.append(row)

    return rows, thresholds


def analyze(rows, thresholds):
    """For each threshold, compute summary stats for cut samples."""
    results = []

    for thresh in thresholds:
        cut_tbl = []

        for row in rows:
            cell = row[thresh].strip()
            tbl = int(row["terminal_branch_length"].strip())

            # Sample is "cut" when its FRS is above the threshold,
            # indicated by the cell starting with '>'
            if cell.startswith(">"):
                cut_tbl.append(tbl)

        n_cut = len(cut_tbl)

        if n_cut == 0:
            mean_tbl   = "N/A"
            median_tbl = "N/A"
            pct_le2    = "N/A"
        else:
            mean_tbl   = round(statistics.mean(cut_tbl), 4)
            median_tbl = round(statistics.median(cut_tbl), 4)
            pct_le2    = round(100 * sum(1 for t in cut_tbl if t <= 2) / n_cut, 2)

        results.append({
            "threshold":              thresh,
            "n_cut":                  n_cut,
            "mean_terminal_bl":       mean_tbl,
            "median_terminal_bl":     median_tbl,
            "pct_cut_with_tbl_le2":   pct_le2,
        })

    return results


def print_results(results, total_samples):
    tsv_header = "\t".join([
        "threshold", "n_cut", "pct_of_total", "mean_TBL", "median_TBL", "pct_cut_TBL_le2"
    ])
    print(tsv_header)

    for r in results:
        pct_total = (
            round(100 * r["n_cut"] / total_samples, 2)
            if total_samples > 0 else "N/A"
        )
        print("\t".join([
            r["threshold"],
            str(r["n_cut"]),
            str(pct_total),
            str(r["mean_terminal_bl"]),
            str(r["median_terminal_bl"]),
            str(r["pct_cut_with_tbl_le2"]),
        ]))


def main():
    filepath = sys.argv[1] if len(sys.argv) > 1 else "metadata.tsv"

    print(f"Reading: {filepath}", file=sys.stderr)
    rows, thresholds = parse_metadata(filepath)
    print(f"Total samples : {len(rows)}", file=sys.stderr)
    print(f"Thresholds    : {thresholds}\n", file=sys.stderr)

    results = analyze(rows, thresholds)
    print_results(results, len(rows))


if __name__ == "__main__":
    main()