To determine the effects of global sampling bias, we will run introduction analysis on a global tree where over represented countries are
 subsampled. 

All commands run from `sampling_bias_analysis`

To subsample overrepresented countries, first we quantify country distribution.
`tail -n +2 ../data/sample_country.tsv | cut -f2 | sort | uniq -c | sort -nr | awk '{print $2"\t"$1}' > data/country_count.tsv`

In `data/country_count.tsv` the 7 countries with the highest counts are:

GBR     14724
CHN     9780
IND     7083
PER     6795
ZAF     6488
AUS     4940
USA     4635

To downsample from these countries, we will make tsv files of samples to be removed. To make a file for removing 10% of the GBR samples we use the
command `python3 scripts/subsample_specific_countries.py --sample_country_tsv ../data/sample_country.tsv --output_file data/prune_samples.tsv --countries GBR --downsample .1`

Running `bash scripts/make_files_prune_trees.sh` will generate pruning files for 10-60% of samples for all 7 countries, and pruning files with 10-60%
of multiple countries. After generating the files, this script will prune the samples from the tree to create downsampled trees. 
 Note that currently downsampling will always take the first XX samples for each country and will not shuffle which samples are cut. Scripts will
have to be modified slightly to generate alternate downsamplings. 

If not using `make_files_prune_trees.sh` to make downsampled trees, the command for downsampling is `matUtils extract -i ../data/mtb.20240912.mask10.pb.gz -s {samples to be pruned} -p -o {path to output tree}`


