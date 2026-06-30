# Identifying a threshold for the pruning of putative coinfections from the final phylogenetic tree

In order to identify a reasonable cutoff for pruning coinfections from the tree, I am simulating coinfected samples and calculating FRS values to identify reasonable values for removal. 

## Step 1 of this process is identifying good samples to pair for coinfection.
 
Get tb metadata to get lineage info `wget https://hgdownload.gi.ucsc.edu/hubs/GCF/000/195/955/GCF_000195955.2/UShER_Mtb_SRA/mtb.20240912.metadata.tsv.gz`. 

**I noticed this metadata contains valencia samples with accessions starting with G instead of proper Biosample accessions. I fixed this with the following commands**
Make a conversion tsv
`join -t $'\t' -1 1 -2 2   <(zcat mtb.20240912.metadata.tsv.gz | cut -f1 | grep G | sort)   <(sort -k2,2 spain_conversion_2col_format_fixed.tsv) > GtoBiosamp_conversion.tsv`
Find and replace Valencia samples with Bioaccessions
`awk -F'\t' 'NR==FNR{map[$1]=$2; next} $1 in map{$1=map[$1]} 1' OFS='\t' GtoBiosamp_conversion.tsv mtb.20240912.metadata.tsv > mtb.20260415.metadata.tsv`

Join FRS and lineage info (easier to do this with old metadata file bc VCFs use Valencia sample ids)
`join -1 1 -2 1 -t $'\t' <(cut -f1,4 data/FRSpercentages_april2026.tsv | sort -k1,1 ) <(cut -f1,10  mtb.20240912.metadata.tsv | sort -k1,1 )  > FRSandlineages.20260415.tsv`
Rename Valencia ids to Biosamples
`awk -F'\t' 'NR==FNR{map[$1]=$2; next} $1 in map{$1=map[$1]} 1' OFS='\t' GtoBiosamp_conversion.tsv FRSandlineages.20260415.tsv > FRSandlineages.20260415_renamed.tsv`


## Identify samples with high FRS values that we can assume are unlikely to be coinfected.

Use `FRSandlineages.20260415_renamed.tsv` to identify good samples to pair for mixing. 
First only use samples with very low FRS proportions 
`awk '$2<=0.03' FRSandlineages.20260415_renamed.tsv  > lowFRSandlineages.tsv`


## Identify high FRS pairs in different lineages so they are different enough to be noticably coinfected. 

Use AI assisted script to randomly pair low FRS samples to other samples cross lineage and within lineages (beyond 2nd order sublineages)
`py scripts/pair_samples.py lowFRSandlineages.tsv pairedformixing.tsv`

## Generate VCFs for samples 

From `pairedformixing.tsv` we randomly selected 100 sample mixes which we store in `sample_pairs_list.txt`

Run pipeline to make and analyze data. `snakemake -s frs_testing.smk --executor slurm --jobs 50 --resources ncbi_api=2 --rerun-incomplete --keep-going`

The pipeline downloads FASTQs for each sample in the final dataset, uses seqtk to mix FASTQs at the appropriate ratio, uses myco to assemble FASTQ and generate VCF, and calculates percentage of low FRS variants in each sample-ratio. 

After pipeline finishes analyze mixed sample stats: 
`python3 scripts/summarize_frs.py /private/groups/corbettlab/lily/frs_testing/frs_matrix.tsv data/summary.tsv`

`data/summary.tsv` summarizes the entirety of the mixed sample FRS experiment. It is essential for future QC decisions.


you can verify samples are mixed with command like this 
`grep "^@" mixed_70_30_R1.fastq | cut -d'.' -f1 | sort | uniq -c`

## After determining mean percentage of low-FRS for each sample mix ratio, look back at the real data.

From the myco pipeline, 135,180 samples passed QC and had data generated for the final tree. The 135,180 sample tree is available at `data/may2026_treebeforepruning_135180samples.pb` From this tree, 137 samples were removed as they were presumed to be *M. cannettii* the pruned samples are at `data/cannettii_samples.tsv` the tree was pruned with the command `matUtils extract -i data/may2026_treebeforepruning_135180samples.pb -s data/cannettii_samples.tsv -p -o data/may2026_treebeforepruning_135180samples.137cannettiiremoved.pb` and the sans-cannettii tree is in `data/may2026_treebeforepruning_135180samples.137cannettiiremoved.pb`. This tree was used for the FRS analysis. 

FRS data was calculated by evaluating the vcfs generated with the following bash script:

`
#!/bin/bash

# Directory containing VCF files
VCF_DIR="" #location of VCFs

echo -e "sample\tvariants\tlow_FRS_variants\tprop_lowFRS_vars"
# Iterate over all .vcf files in the directory
for inputfile in "$VCF_DIR"/*.vcf; do
	# Extract the filename without the directory and extension
	filename=$(basename "$inputfile" .vcf)
	#get total number of variants in vcf
	tot=$(cat "$inputfile" | grep -v "#" | wc -l)
	# Run bcftools filter command
	# identify number of variants below FRS threshold
	filt=$(bcftools filter -i 'FRS<.85' "$inputfile" | grep -v '#' | wc -l)
	#calculate percentage of variants below FRS threshold
	result=$(echo "scale=4; $filt/$tot" | bc)
	echo -e "${filename}\t${tot}\t${filt}\t${result}"

	# Optional: Print status
	#echo "Filtered $inputfile -> $outputfile"
done
`

It detects the number of "low-FRS" (<0.85) variants in a vcf compared to the total number of variants. **NOTE: due to storage demands VCF files are not publicly available. Pls contact the corresponding author for further information.**

The results of this script are found in `data/FRSpercents_unpruned_may2026_FINAL.tsv` this file contains an FRS score for all 135,180 samples on the tree (cannettii included).

The FRS values for samples minus cannettii can be found here `data/FRSpercents_unpruned_may2026_FINAL_nocannettii.tsv`. 

In order to visualize samples with each FRS value, we need to make metadata. 
`python3 scripts/make_metadata_file.py -s data/FRSpercents_unpruned_may2026_FINAL_nocannettii.tsv -t data/summary.tsv -o data/FRSmetadata_nocannettii_FINAL.tsv -pb data/may2026_treebeforepruning_135180samples.canettiremoved.pb`

We then make the visual with 
`usher_to_taxonium -i data/may2026_treebeforepruning_135180samples.137cannettiiremoved.pb -m data/FRSmetadata_nocannettii_FINAL.tsv -c strain,frs,0.4,0.39,0.38,0.37,0.36,0.34,0.32,0.31,0.28,0.26,0.24,0.22,0.19,0.17,0.15,0.12,0.09,0.07,0.05,0.04,terminal_branch_length,root_distance -o data/may2026_treebeforepruning_135180samples.137cannettiiremoved.annotated.jsonl.gz`

`data/may2026_treebeforepruning_135180samples.137cannettiiremoved.annotated.jsonl.gz` can be input to taxonium.org to visualize the cutoffs. 

`python scripts/metadata_analysis.py data/FRSmetadata_nocannettii_FINAL.tsv > data/cutoff_statistics.tsv` will identify the number of samples removed by each threshold and identify the average and median terminal branch length as well as the number of samples with a terminal branch length of 2 or less. 

To get a list of samples to prune from the tree based on FRS value > 0.19 `awk '$2 > 0.19' data/FRSmetadata_nocannettii_FINAL.tsv > data/FRSgreaterthan.19.tsv`

prune these samples from the data `matUtils extract -i data/may2026_treebeforepruning_135180samples.137cannettiiremoved.pb -s data/FRSgreaterthan.19.tsv -p -o data/may2026_treebeforepruning_135180samples.137cannettiiremoved.5631samplesFRSpruned.pb`

After you decide that 0.19 is the best cutoff for FRS you need to make one more visual.
`python3 scripts/annotate_manual_prunes.py --prune-list  ../Supplemental/SupplementalTable3.tsv  --large-file  data/FRSmetadata_nocannettii_FINAL.tsv --output data/FRSmetadata_nocannettii_FINAL_plusmanual.tsv`

`usher_to_taxonium -i data/may2026_treebeforepruning_135180samples.137cannettiiremoved.pb -m data/FRSmetadata_nocannettii_FINAL_plusmanual.tsv -c strain,frs,0.4,0.39,0.38,0.37,0.36,0.34,0.32,0.31,0.28,0.26,0.24,0.22,0.19,0.17,0.15,0.12,0.09,0.07,0.05,0.04,terminal_branch_length,root_distance,filter -o data/may2026_treebeforepruning_135180samples.137cannettiiremoved.fullannotation.jsonl.gz`







