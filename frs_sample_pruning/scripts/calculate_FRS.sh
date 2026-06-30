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