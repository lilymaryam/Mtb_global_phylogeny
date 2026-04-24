This analysis was used to generate code for Supplemental Figure 2 (maybe an additional figure addressing the impact of unequal global sampling)

From `introduction_analysis` make a data dir `mkdir data`

Download these files into `introduction_analysis/data`

Tree: https://hgdownload.gi.ucsc.edu/hubs/GCF/000/195/955/GCF_000195955.2/UShER_Mtb_SRA/
 
`mtb.20240912.mask10.pb.gz              2025-07-07 16:23   17M`

metadata: https://hgdownload.gi.ucsc.edu/hubs/GCF/000/195/955/GCF_000195955.2/UShER_Mtb_SRA/

`mtb.20240912.metadata.tsv.gz           2026-04-16 11:45  1.8M`

Run all commands from `introduction_analysis` directory 

extract relevant metadata (sample_index    country region  tbprofiler_lineage_usher)
`zcat data/mtb.20240912.metadata.tsv.gz | cut -f1,3,4,10 > data/relevant_metadata.tsv`

*dont really need main lineage but still fine*
get sample country and main lineage
`cut -f1,2,4 data/relevant_metadata.tsv | perl -pe 's/lineage/L/g' | awk -F'\t' 'BEGIN {OFS="\t"} {gsub(/\..*/, "", $3); print $1, $2, $3}' > data/sample_country_lineage.tsv`

get sample and country
`cut -f1,2 data/sample_country_lineage.tsv > data/sample_country.tsv`

introduction analysis
`matUtils introduce -i data/mtb.20240912.mask10.pb.gz -s data/sample_country.tsv -o data/introduce_output.tsv`

generate figure/breakdown of main introductions
`python3 scripts/analyze_intros.py`

Will generate `Full_Tree_Introduction_Analysis.png` in `introduction_analysis` This is Supplemental Figure 2. 

