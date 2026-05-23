# Improving and automating the curation of ancestrally placed samples on the tree
In order to prove my method for pruning coinfections is valid I am simulating coinfected samples and determining valid cutoffs for removal. 
Step 1 of this process is identifying good samples to pair for coinfection. This is a 2 part process.
Part 1: 
Identify samples with high FRS values that we can assume are unlikely to be coinfected. 
Part 2:
Identify high FRS pairs in different lineages so they are different enough to be noticably coinfected. 

Step 1:
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

Use `FRSandlineages.20260415_renamed.tsv` to identify good samples to pair for mixing. 
First only use samples with very low FRS proportions 
`awk '$2<=0.03' FRSandlineages.20260415_renamed.tsv  > lowFRSandlineages.tsv`

Use AI assisted script to randomly pair low FRS samples to other samples cross lineage and within lineages (beyond 2nd order sublineages)
`py scripts/pair_samples.py lowFRSandlineages.tsv pairedformixing.tsv`

From `pairedformixing.tsv` we randomly selected 100 sample mixes which we store in `sample_pairs_list.txt`

Run pipeline to make and analyze data. `snakemake -s frs_testing.smk --executor slurm --jobs 50 --resources ncbi_api=2 --rerun-incomplete --keep-going`

Analyze mixed sample stats 
`python3 scripts/summarize_frs.py /private/groups/corbettlab/lily/frs_testing/frs_matrix.tsv summary`

verify samples are mixed with command like this 
`grep "^@" mixed_70_30_R1.fastq | cut -d'.' -f1 | sort | uniq -c`


## mixed FASTQs are run through myco pipeline
`wget https://raw.githubusercontent.com/aofarrel/myco/7.0.6/myco_raw.wdl`
create template for inputs to myco_raw
`miniwdl input-template myco_raw.wdl > inputs_template.json`
