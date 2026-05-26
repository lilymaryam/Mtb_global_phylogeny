import argparse
import bte
import subprocess

#working on this 
'''
def make_jsonl(out, pb, thresholds):
    cmd = f"usher_to_taxonium -i data/originaltree.beforecoinfprune.aftercanprune.20250319.rootadded.reroot.col.renamed.pb -m test -c strain,0.4,0.39,0.38,0.37,0.36,0.34,0.32,0.31,0.28,0.26,0.24,0.22,0.19,0.17,0.15,0.12,0.09,0.07,0.05,0.04 -o test.jsonl.gz"
'''
    
'''
def get_path_length(samples, pb_path, path_lengths):
    # Use bte to load the tree and get the path length
    tree = bte.MATree(pb_path)
    # Assuming we want the path length from root to each leaf
    with open(samples) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 4:
                sample_id = parts[0]
                frs = float(parts[3])
                if tree.get_node(sample_id):
                    print(tree.get_node(sample_id))
                else:
                    print(f"Sample {sample_id} not found in the tree") # Ensure the sample is in the tree
                #path_lengths[sample_id] = path_length
    #for leaf in tree.leaves():
    #    path_length = tree.get_path_length(leaf)
    #    path_lengths[leaf.name] = path_length
    #return path_lengths
'''

def get_path_length(samples, pb_path, path_lengths):
    # Use bte to load the tree and get the path length
    tree = bte.MATree(pb_path)
    
    with open(samples) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 4:
                sample_id = parts[0]
                frs = float(parts[3])
                
                try:
                    node = tree.get_node(sample_id)
                    # Get path length from root to this node
                    path_length = 0
                    current = node
                    while current.parent:
                        path_length += current.branch_length if current.branch_length else 0
                        current = current.parent
                    path_lengths[sample_id] = path_length
                except ValueError:
                    # Node not found in tree
                    print(f"Sample {sample_id} not found in the tree")
    
    return path_lengths

def get_thresholds(threshold_file):
    thresholds = []
    with open(threshold_file) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                #ratio = parts[0]
                #mean_frs = float(parts[1]).round(2)
                mean_frs = float(f"{float(parts[1]):.2f}")
                #thresholds.append((ratio, mean_frs))
                thresholds.append(mean_frs)
    return thresholds

def new_metadata(samples_file, thresholds, output, pb_path):
    tree = bte.MATree(pb_path)

    with open(output, 'w') as out:
        above_thresholds = [str(th) for th in thresholds]
        header = '\t'.join(['strain', 'frs'] + above_thresholds + ['terminal_branch_length'])
        out.write(f"{header}\n")
        with open(samples_file) as f:
            #header = f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 4:
                    sample_id = parts[0]
                    frs = float(parts[3])
                    thresh = []
                    path_length = "NA"
                    actual_threshold = None
                    try:
                        node = tree.get_node(sample_id)
                        # Get path length from root to this node
                        path_length = len(node.mutations)
                        #i dont actually want path just branch lenght 
                        #path_length = 0
                        #current = node
                        #while current.parent:
                        #    path_length += current.branch_length if current.branch_length else 0
                        #    current = current.parent
                        
                        #path_lengths[sample_id] = path_length
                    except ValueError:
                        # Node not found in tree
                        print(f"Sample {sample_id} not found in the tree")
                    for th in thresholds:
                        if frs > th:
                            thresh.append(f">{th:.2f}")
                            actual_threshold = th
                            #print(f"{sample_id}\t{frs:.4f}\t\t{th:.2f}")
                        else:
                            thresh.append(f"<={th:.2f}")
                            #print(f"{sample_id}\t{frs:.4f}\t\t{th:.2f}")
                    out.write(f"{sample_id}\t{frs:.4f}\t" + '\t'.join(thresh) + f"\t{path_length}" + "\n")
                    #lineage = parts[2]
                    # Determine which thresholds this sample is above
                    #above_thresholds = [str(th) for th in thresholds if frs > th]
                    #print(f"{sample_id}\t{frs:.4f}\t\t{','.join(above_thresholds)}")

def main():
    parser = argparse.ArgumentParser(description='Count samples with FRS above each ratio threshold')
    parser.add_argument('-s', '--samples', required=True, help='Sample file (sample_id, frs, lineage)')
    parser.add_argument('-t', '--thresholds', required=True, help='Summary TSV from FRS analysis (ratio, mean, stddev, n)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('-pb', '--pb_path', required=True, help='Output TSV file')
    #parser.add_argument('-j', '--jsonl_path', required=True, help='Output JSONL file')

    args = parser.parse_args()

    threshold_nums = get_thresholds(args.thresholds)
    new_metadata(args.samples, threshold_nums, args.output, args.pb_path)
    #paths = {}
    #get_path_length(args.samples, args.pb_path, paths)


if __name__ == "__main__":
    main()

