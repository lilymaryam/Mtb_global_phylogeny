# This script obtains a pairwise-distance matrix for a given phylogeny

library(castor)
library(ape)

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)

# Read tree in newick format and get a list of the name of the samples
tree_usher <- read.tree(args[1])
samples <- tree_usher$tip.label

print("CALCULATING PAIRWISE DISTANCES")
distances_from_tree <- get_all_pairwise_distances(tree_usher, only_clades = samples)
colnames(distances_from_tree) <- samples
rownames(distances_from_tree) <- samples

print("SAVING DISTANCES IN A MATRIX")
write.csv(distances_from_tree, paste("distancematrix_",tree_usher, sep=""), row.names = T)
