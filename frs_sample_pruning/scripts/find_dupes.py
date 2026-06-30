import pandas as pd
import sys

def find_duplicate_samples(samples_file, conversion_file, output_file=None):
    """
    Find G samples that are duplicates (already present as SAMEA names).
    
    Args:
        samples_file: Path to file with sample, parsimony, parent_id columns
        conversion_file: Path to file with G to SAMEA name mappings
        output_file: Optional path to save the list of G samples to remove
    """
    
    # Read the samples file
    samples_df = pd.read_csv(samples_file, sep='\s+')
    
    # Read the conversion mapping file
    conversion_df = pd.read_csv(conversion_file, sep='\s+', header=None, 
                                names=['G_name', 'SAMEA_name'])
    
    # Create a mapping dictionary
    name_mapping = dict(zip(conversion_df['G_name'], conversion_df['SAMEA_name']))
    
    # Get all samples from the samples list
    all_samples = set(samples_df['sample'])
    
    # Find G samples that should be removed (both G and SAMEA versions exist)
    samples_to_remove = []
    for g_name, samea_name in name_mapping.items():
        if g_name in all_samples and samea_name in all_samples:
            samples_to_remove.append(g_name)
    
    # Print results
    if samples_to_remove:
        print(f"Found {len(samples_to_remove)} G samples to remove (duplicates):\n")
        for sample in samples_to_remove:
            print(sample)
        
        # Save to file if requested
        if output_file:
            with open(output_file, 'w') as f:
                f.write('\n'.join(samples_to_remove) + '\n')
            print(f"\nList saved to: {output_file}")
    else:
        print("No duplicate G samples found")
    
    return samples_to_remove

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <samples_file> <conversion_file> [output_file]")
        print("\nExample:")
        print("  python script.py samples.txt conversions.txt")
        print("  python script.py samples.txt conversions.txt samples_to_remove.txt")
        sys.exit(1)
    
    samples_file = sys.argv[1]
    conversion_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    find_duplicate_samples(samples_file, conversion_file, output_file)