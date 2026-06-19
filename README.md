# This repository contains scripts for analysis and visualization of data for the manuscript UShER-TB. 
### Transmission cluster assignment in local and global datasets
**cluster_benchmarking_metrics.R**: R script used to calculate adjusted Rand-index, precision and recall to compare clusters obtained with the gold-standard approach and with usher-sampled placement. 

**placement_comparison_and_benchmarking_files**: zip file containing all the files used to perform the initial benchmarking and the following cluster comparisons 

**compare_clusters.R**: R script used to compare clusters identified by standard approach for cluster identification (https://gitlab.com/tbgenomicsunit/ThePipeline) using a 10 SNPs threshold and clusters identified by usher-sampled placement.

**distance_matrix_from_phylogeny.R**: R script used to obtain the pairwise genetic distance between all samples. This matrix is input for distance_from_R_matrix.py

**distance_from_R_matrix.py**: Python script that obtains the clusters for a given SNP threshold (10 in our case) which it uses for comparison against the input pairwise genetic distance matrix from distance_matrix_from_phylogeny.R.

**cohens_kappa.R**: R script to calculate the Cohen's Kappa index to evaluate usher-sampled accuracy of assigning newly-placed samples to their corresponding cluster 

### Ancient MTBC genomes placement with usher-sampled 
**vcf_to_diff_script.py**: Python script converts vcf to MAPLE format. This version requires a bed file of commonly masked regions of the MTB genome found at https://github.com/iqbal-lab-org/cryptic_tb_callable_mask/blob/master/R00000039_repregions.bed and provides an option to provide a bed file of low coverage regions in the sample genome for masking. 

### Evaluation of the accuracy of UShER for public health inquiry
**nrr_analysis_pipeline/**: this directory contains a snakemake pipeline that compares the accuracy of the phylogeny built by usher-sampled with simulated data from Alisim (nrr_analysis_pipeline/data/truth.30000.col.pb)(see Methods) to the truth phylogeny (nrr_analysis_pipeline/data/truth.30000.col.pb). To run this pipeline, naviagate to the dir and run `snakmake -s nrr_analysis --cores {number of cores}`. This will result with nrr_visual.svg, which is Extended Data Figure 1 in the paper. 

**rfdendro.py**: calculates the normalized robinson foulds distance between 2 bifurcating newick trees WITH IDENTICAL TAXA
To calculate the normalized rf distance between the simulated and truth trees use `data/*.30000.nwk`

### _De novo_ phylogeny reconstruction with UShER-TB
Once the Variant Calling has been performed using the _myco_ pipeline, the tree reconstruction consists of an easy, 2-step pipeline:

1. First, we need to prepare the base tree with the reference by generating a pb (protobuf) file, which will act as a starting point for sample placement:

```echo "();" > base_tree.nwk```

```echo ">ref" > ref.diff"```

```usher-sampled -t base_tree.nwk --diff ref.diff --ref <reference.fasta> -o ref.pb```

2. Then, we add the rest of the samples to the ref.pb tree, generating a new, complete pb tree. To generate the newick file, we add the `-u` parameter:

```usher-sampled -i <ref.pb> --diff <concatenated_diff_files> --ref <reference.fasta> -o <output_tree.pb> -u -d <output_directory>```


### Sample placement in an existing tree using UShER-TB
The basic command needed to add new sequences to an existing phylogeny is as follows:

`usher-sampled -i <tree_in_pb_format> -o <output_name.pb> --diff <diff_file> -d <output_directory>`

Where:
- The tree in pb (protobuf) format can be obtained by reconstructing a phylogeny _de novo_ or annotating an existing newick file with the mutations in VCF format with UShER-TB.
- The diff files are obtained for each sample individually using the _myco_ pipeline. Several files can be concatenated in a single diff file for use as input for UShER-TB.









