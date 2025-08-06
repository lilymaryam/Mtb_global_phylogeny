import dendropy
from dendropy.calculate import treecompare
import argparse

def normalized_rf_distance(tree_file1, tree_file2):
    # Load trees
  shared_taxon_namespace = dendropy.TaxonNamespace()
  # Load trees with the shared TaxonNamespace
  tree1 = dendropy.Tree.get(path=tree_file1, schema="newick", taxon_namespace=shared_taxon_namespace)
  tree2 = dendropy.Tree.get(path=tree_file2, schema="newick", taxon_namespace=shared_taxon_namespace)

  # Calculate unweighted RF distance
  unweighted_rf = treecompare.unweighted_robinson_foulds_distance(tree1, tree2)

  # Count bipartitions
  tree1.encode_bipartitions()
  tree2.encode_bipartitions()
  bipartitions_tree1 = len(tree1.bipartition_encoding)
  bipartitions_tree2 = len(tree2.bipartition_encoding)
  #print(tree1.is_bipartitions_encoded)  # Should return True
  #print(tree2.is_bipartitions_encoded)
  # Calculate maximum RF distance
  #for bipartition in tree1.bipartition_encoding:
  #  print(bipartition)
  internal_edges = len(tree1.internal_edges())
  print(internal_edges)
  max_rf = bipartitions_tree1 + bipartitions_tree2

  # Normalize RF distance
  normalized_rf = unweighted_rf / max_rf if max_rf > 0 else 0
  
  return unweighted_rf, max_rf, normalized_rf

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate normalized Robinson-Foulds distance between two phylogenetic trees.")
    parser.add_argument("-t1", "--tree1", required=True, help="Path to the first tree file (Newick format)")
    parser.add_argument("-t2", "--tree2", required=True, help="Path to the second tree file (Newick format)")
    return parser.parse_args()

# Example usage
def main():
    args = parse_args()
    
    unweighted_rf, max_rf, normalized_rf = normalized_rf_distance(args.tree1, args.tree2)
    
    print(f"Unweighted RF Distance: {unweighted_rf}")
    print(f"Maximum RF Distance: {max_rf}")
    print(f"Normalized RF Distance: {normalized_rf}")

if __name__ == "__main__":
    main()