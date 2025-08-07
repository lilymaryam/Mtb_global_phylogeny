# This repository contains scripts for analysis and visualization of data for the manuscript UShER-TB. 
### Transmission cluster assignment in local and global datasets
compare_clusters.R: R script used to compare clusters identified by standard approach for cluster identification (https://gitlab.com/tbgenomicsunit/ThePipeline) using a 10 SNPs threshold and clusters identified by usher-sampled placement.

distance_matrix_from_phylogeny.R: R script used to obtain the pairwise genetic distance between all samples. This matrix is input for distance_from_R_matrix.py

distance_from_R_matrix.py: Python script that obtains the clusters for a given SNP threshold (10 in our case) which it uses for comparison against the input pairwise genetic distance matrix from distance_matrix_from_phylogeny.R.

### Ancient MTBC genomes placement with usher-sampled 
vcf_to_diff_script.py: Python script converts vcf to MAPLE format. This version requires a bed file of commonly masked regions of the MTB genome found at https://github.com/iqbal-lab-org/cryptic_tb_callable_mask/blob/master/R00000039_repregions.bed and provides an option to provide a bed file of low coverage regions in the sample genome for masking. 

### Evaluation of the accuracy of UShER for public health inquiry
nrr_analysis_pipeline/: this directory contains a snakemake pipeline that compares the accuracy of the phylogeny built by usher-sampled with simulated data from Alisim (nrr_analysis_pipeline/data/truth.30000.col.pb)(see Methods) to the truth phylogeny (nrr_analysis_pipeline/data/truth.30000.col.pb). To run this pipeline, naviagate to the dir and run `snakmake -s nrr_analysis --cores {number of cores}`. This will result with nrr_visual.svg, which is Extended Data Figure 1 in the paper. 

rfdendro.py: calculates the normalized robinson foulds distance between 2 bifurcating newick trees WITH IDENTICAL TAXA
To calculate the normalized rf distance between the simulated and truth trees use `data/*.30000.nwk`










