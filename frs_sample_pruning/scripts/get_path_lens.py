#!/usr/bin/env python3
import sys
from collections import defaultdict
import statistics

def extract_cutoff_value(cutoff_str):
    """Extract numeric value from cutoff string like '<=0.40' or '>0.05'"""
    return float(cutoff_str.replace('<=', '').replace('>', ''))

def find_best_cutoff(cutoffs):
    """Find the first '>' cutoff (or None if all are '<=')"""
    for cutoff in cutoffs:
        if cutoff.startswith('>'):
            return cutoff
    # If all are '<=', return None (no cutoff)
    return None

def process_data(input_file):
    """Process data and group last column by best cutoff"""
    
    cutoff_data = defaultdict(list)
    skipped_na = 0
    skipped_no_cutoff = 0
    
    with open(input_file, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split()
            
            if len(parts) < 23:
                continue
            
            sample_id = parts[0]
            real_value = float(parts[1])
            cutoffs = parts[2:22]  # 20 cutoff columns
            
            # Skip if last value is NA
            if parts[22] == 'NA':
                skipped_na += 1
                continue
            
            last_value = float(parts[22])
            
            # Find best cutoff
            best_cutoff = find_best_cutoff(cutoffs)
            
            # Skip if no cutoff found (all are '<=')
            if best_cutoff is None:
                skipped_no_cutoff += 1
                continue
            
            # Append last column value to this cutoff's list
            cutoff_data[best_cutoff].append(last_value)
            
            print(f"{sample_id}: real={real_value:.4f}, best_cutoff={best_cutoff}, value={last_value}")
    
    if skipped_na > 0:
        print(f"\nSkipped {skipped_na} samples with NA values")
    if skipped_no_cutoff > 0:
        print(f"Skipped {skipped_no_cutoff} samples with no cutoff (all <=)")
    
    return cutoff_data

'''
def write_results(cutoff_data):
    """Write summary of results to file"""
    with open("cutoff_results.txt", 'w') as out:
        out.write("\n" + "="*60)
        out.write("\nRESULTS: Values grouped by best cutoff")
        out.write("\n" + "="*60)
        
        # Sort by cutoff value for readable output
        for cutoff in sorted(cutoff_data.keys(), key=extract_cutoff_value, reverse=True):
            values = cutoff_data[cutoff]
            out.write(f"\n{cutoff}: {len(values)} samples")
            out.write(f"\n  Values: {values[:10]}{'...' if len(values) > 10 else ''}")
            out.write(f"\n  Mean: {sum(values)/len(values):.2f}")
            out.write(f"\n  Median: {statistics.median(values):.2f}")
            out.write(f"\n  Min: {min(values):.2f}, Max: {max(values):.2f}")
            '''

def write_results(cutoff_data):
    """Write summary of results to TSV file"""
    with open("cutoff_results.tsv", 'w') as out:
        # Write header
        out.write("cutoff\tnum_samples\tmean\tmedian\tmin\tmax\tpct_pathlen_le2\n")
        
        # Sort by cutoff value for readable output
        for cutoff in sorted(cutoff_data.keys(), key=extract_cutoff_value, reverse=True):
            values = cutoff_data[cutoff]
            pct_le2 = 100 * sum(1 for v in values if v <= 5) / len(values)
            out.write(f"{cutoff}\t{len(values)}\t{sum(values)/len(values):.2f}\t{statistics.median(values):.2f}\t{min(values):.2f}\t{max(values):.2f}\t{pct_le2:.1f}\n")

'''
def write_results(cutoff_data):
    """Write summary of results to TSV file"""
    with open("cutoff_results.tsv", 'w') as out:
        # Write header
        out.write("cutoff\tnum_samples\tmean\tmedian\tmin\tmax\n")
        
        # Sort by cutoff value for readable output
        for cutoff in sorted(cutoff_data.keys(), key=extract_cutoff_value, reverse=True):
            values = cutoff_data[cutoff]
            out.write(f"{cutoff}\t{len(values)}\t{sum(values)/len(values):.2f}\t{statistics.median(values):.2f}\t{min(values):.2f}\t{max(values):.2f}\n")

    
    # Sort by cutoff value for readable output
    for cutoff in sorted(cutoff_data.keys(), key=extract_cutoff_value, reverse=True):
        values = cutoff_data[cutoff]
        print(f"\n{cutoff}: {len(values)} samples")
        print(f"  Values: {values[:10]}{'...' if len(values) > 10 else ''}")
        print(f"  Mean: {sum(values)/len(values):.2f}")
        print(f"  Median: {statistics.median(values):.2f}")
        print(f"  Min: {min(values):.2f}, Max: {max(values):.2f}")
        '''

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    cutoff_data = process_data(input_file)
    write_results(cutoff_data)
    
    # Save to file
    '''
    with open("cutoff_results.txt", 'w') as out:
        for cutoff in sorted(cutoff_data.keys(), key=extract_cutoff_value, reverse=True):
            values = cutoff_data[cutoff]
            out.write(f"{cutoff}\t{','.join(map(str, values))}\n")
            '''
    
    print("\nResults saved to cutoff_results.txt")