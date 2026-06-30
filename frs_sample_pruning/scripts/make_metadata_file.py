import argparse
import bte
import subprocess

def get_thresholds(threshold_file):
    thresholds = []
    with open(threshold_file) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                mean_frs = float(f"{float(parts[1]):.2f}")
                thresholds.append(mean_frs)
    return thresholds

def new_metadata(samples_file, thresholds, output, pb_path):
    tree = bte.MATree(pb_path)

    with open(output, 'w') as out:
        above_thresholds = [str(th) for th in thresholds]
        header = '\t'.join(['strain', 'frs'] + above_thresholds + ['terminal_branch_length', 'root_distance'])
        out.write(f"{header}\n")
        with open(samples_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 4:
                    sample_id = parts[0]
                    frs = float(parts[3])
                    thresh = []
                    branch_length = "NA"
                    root_distance = "NA"

                    try:
                        node = tree.get_node(sample_id)
                        branch_length = len(node.mutations)
                        root_distance = 0
                        current = node
                        while current.parent:
                            root_distance += len(current.mutations)
                            current = current.parent
                    except ValueError:
                        print(f"Sample {sample_id} not found in the tree")

                    for th in thresholds:
                        if frs > th:
                            thresh.append(f">{th:.2f}")
                        else:
                            thresh.append(f"<={th:.2f}")

                    out.write(f"{sample_id}\t{frs:.4f}\t" + '\t'.join(thresh) + f"\t{branch_length}\t{root_distance}\n")

def main():
    parser = argparse.ArgumentParser(description='Count samples with FRS above each ratio threshold')
    parser.add_argument('-s', '--samples', required=True, help='Sample file (sample_id, frs, lineage)')
    parser.add_argument('-t', '--thresholds', required=True, help='Summary TSV from FRS analysis (ratio, mean, stddev, n)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('-pb', '--pb_path', required=True, help='Path to .pb tree file')

    args = parser.parse_args()

    threshold_nums = get_thresholds(args.thresholds)
    new_metadata(args.samples, threshold_nums, args.output, args.pb_path)

if __name__ == "__main__":
    main()
