#!/bin/bash

# Define arrays
countries=(GBR CHN IND PER ZAF AUS USA)
downsamples=(0.1 0.2 0.3 0.4 0.5 0.6)

# Check we're in the right directory
if [[ $(basename "$PWD") != "sampling_bias_analysis" ]]; then
    echo "Error: This script must be run from the sampling_bias_analysis directory"
    echo "Current directory: $(basename "$PWD")"
    exit 1
fi

if [[ ! -f "scripts/subsample_specific_countries.py" ]]; then
    echo "Error: Cannot find scripts/subsample_specific_countries.py"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p data

echo "Running individual country subsampling..."

# Loop through each country separately
for country in "${countries[@]}"; do
    for downsample in "${downsamples[@]}"; do
        # Convert downsample to integer for filename (e.g., 0.1 -> 10)
        pct=$(echo "$downsample * 100" | bc | cut -d. -f1)
        
        samples_file="data/prune_samples_${country}_${pct}pct.tsv"
        tree_file="data/tree_${country}_${pct}pct.pb"
        
        echo "Processing $country at ${pct}% downsample..."
        python3 scripts/subsample_specific_countries.py \
            --sample_country_tsv ../data/sample_country.tsv \
            --output_file "$samples_file" \
            --countries "$country" \
            --downsample "$downsample"
        
        echo "Extracting tree for $country at ${pct}%..."
        matUtils extract \
            -i ../data/mtb.20240912.mask10.pb.gz \
            -s "$samples_file" \
            -p \
            -o "$tree_file"
    done
done

echo "Running combined country subsampling..."

# Loop through downsamples with all countries together
for downsample in "${downsamples[@]}"; do
    pct=$(echo "$downsample * 100" | bc | cut -d. -f1)
    
    samples_file="data/prune_samples_ALL_${pct}pct.tsv"
    tree_file="data/tree_ALL_${pct}pct.pb"
    
    # Join countries with commas
    countries_str=$(IFS=,; echo "${countries[*]}")
    
    echo "Processing all countries at ${pct}% downsample..."
    python3 scripts/subsample_specific_countries.py \
        --sample_country_tsv ../data/sample_country.tsv \
        --output_file "$samples_file" \
        --countries "$countries_str" \
        --downsample "$downsample"
    
    echo "Extracting tree for all countries at ${pct}%..."
    matUtils extract \
        -i ../data/mtb.20240912.mask10.pb.gz \
        -s "$samples_file" \
        -p \
        -o "$tree_file"
done

echo "Done!"