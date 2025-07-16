# This script converts the distance matrix obtained using the distance_matrix_from_phylogeny.R script into a genetic distances file (format Sample\tDistance\tSample2) filtering by a maximum SNP distance of 12. 

import pandas as pd
import numpy as np
import sys

## First, we import the distance matrix obtained with R script distance_matrix_from_phylogeny.R 

print("Importing distance matrix")
df = pd.read_csv(sys.argv[1], index_col=0)

print("Extracting upper matrix")
upper_tri_indices = np.triu_indices_from(df, k=1)

upper_df = pd.DataFrame({
    'Strain1': df.index[upper_tri_indices[0]],
    'Strain2': df.columns[upper_tri_indices[1]],
    'Diff': df.values[upper_tri_indices]
})

# Then we convert the matrix to long-format for efficiency
print("Converting matrix to long format")

long_df = df.stack().reset_index()
long_df.columns = ['Strain1', 'Strain2', 'Diff']
long_df['Diff'] = long_df['Diff'].astype(int)

print("Filtering by 12 SNPs")
long_df = long_df[long_df['Diff'] <= 12]
print(long_df.head())

print("Filtering self-comparisons")
# Filter out self-comparisons if necessary (optional)
long_df = long_df[long_df['Strain1'] != long_df['Strain2']]
print(long_df.head())

print("Filtering repeated combinations")
print("Step 1")
# Sort Sample1 and Sample2 for each row to handle undirected pairs
long_df['SortedSamples'] = long_df.apply(lambda x: tuple(sorted([x['Strain1'], x['Strain2']])), axis=1)
print(long_df.head())

print("Step 2")
# Drop duplicate pairs based on SortedSamples
long_df = long_df.drop_duplicates(subset=['SortedSamples'])

# Remove the SortedSamples helper column
long_df = long_df.drop(columns=['SortedSamples'])

print("Filtering by 12 SNPs and ordering")
new_order = ['Strain1', 'Diff', 'Strain2']
long_df = long_df[new_order]
print(long_df.head())

print("Sorting and saving")
sorted_filtered_df = long_df.sort_values(by='Diff')
print(long_df.head())
sorted_filtered_df.to_csv('{0}_genetic_distances_12snps.tsv'.format(sys.argv[1]), sep='\t', index=False)
print("DONE")
