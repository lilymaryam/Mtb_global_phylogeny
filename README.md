# This repository contains scripts for analysis and visualization of data for the manuscript UShER-TB. 
## Results

## Methods
### Transmission cluster assignment in local and global datasets
compare_clusters.R: R script used to compare clusters identified by standard approach for cluster identification (https://gitlab.com/tbgenomicsunit/ThePipeline) using a 10 SNPs threshold and clusters identified by usher-sampled placement.

distance_matrix_from_phylogeny.R: R script used to obtain the pairwise genetic distance between all samples. This matrix is input for distance_from_R_matrix.py

distance_from_R_matrix.py: Python script that obtains the clusters for a given SNP threshold (10 in our case) which it uses for comparison against the input pairwise genetic distance matrix from distance_matrix_from_phylogeny.R.

### Ancient MTBC genomes placement with usher-sampled 
vcf_to_diff_script.py: Python script converts vcf to MAPLE format. This version requires a bed file of commonly masked regions of the MTB genome found at https://github.com/iqbal-lab-org/cryptic_tb_callable_mask/blob/master/R00000039_repregions.bed and provides an option to provide a bed file of low coverage regions in the sample genome for masking. 


