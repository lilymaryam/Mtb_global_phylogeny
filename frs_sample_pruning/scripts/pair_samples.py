#!/usr/bin/env python3

import sys
import random
from collections import defaultdict

def read_samples(filename):
    """Read samples from TSV file and group by lineage"""
    lineage_groups = defaultdict(list)
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split('\t')
            if len(parts) >= 3:
                sample_id = parts[0]
                mixed_value = float(parts[1])
                lineage = parts[2]
                
                lineage_groups[lineage].append({
                    'id': sample_id,
                    'mixed_value': mixed_value,
                    'lineage': lineage
                })
    
    return lineage_groups

def get_lineage_hierarchy(lineage):
    """Extract lineage hierarchy components"""
    # Handle different lineage formats
    if lineage.startswith('lineage'):
        # lineage2.2.1 -> ['lineage2', 'lineage2.2', 'lineage2.2.1']
        parts = lineage.split('.')
        hierarchy = []
        current = parts[0]
        hierarchy.append(current)
        for part in parts[1:]:
            current += '.' + part
            hierarchy.append(current)
        return hierarchy
    elif lineage.startswith('La'):
        # La1.8.1 -> ['La1', 'La1.8', 'La1.8.1']
        parts = lineage.split('.')
        hierarchy = []
        current = parts[0]
        hierarchy.append(current)
        for part in parts[1:]:
            current += '.' + part
            hierarchy.append(current)
        return hierarchy
    else:
        # Single level lineage
        return [lineage]

def can_pair_lineages(lineage1, lineage2):
    """Check if two lineages can be paired based on hierarchy rules"""
    hier1 = get_lineage_hierarchy(lineage1)
    hier2 = get_lineage_hierarchy(lineage2)
    
    # Get major lineage (first level)
    major1 = hier1[0] if hier1 else lineage1
    major2 = hier2[0] if hier2 else lineage2
    
    # Get second-order sublineage (second level if exists)
    second_order1 = hier1[1] if len(hier1) > 1 else hier1[0] if hier1 else lineage1
    second_order2 = hier2[1] if len(hier2) > 1 else hier2[0] if hier2 else lineage2
    
    # Rule 1: Cannot share same second-order sublineage
    if second_order1 == second_order2:
        return False, "same_second_order"
    
    # Rule 2: Prefer cross-major-lineage pairs
    cross_major = (major1 != major2)
    
    return True, "cross_major" if cross_major else "within_major"

def pair_samples(lineage_groups, max_pairs_per_lineage=None):
    """Create pairs of samples from different lineages with balanced distribution"""
    lineages = list(lineage_groups.keys())
    pairs = []
    used_samples = set()
    
    print(f"Found lineages: {lineages}")
    print(f"Sample counts per lineage:")
    for lineage in lineages:
        print(f"  {lineage}: {len(lineage_groups[lineage])}")
    print()
    
    # Get all valid lineage combinations
    cross_major_combos = []
    within_major_combos = []
    
    for i in range(len(lineages)):
        for j in range(i + 1, len(lineages)):
            lineage1 = lineages[i]
            lineage2 = lineages[j]
            can_pair, pair_type = can_pair_lineages(lineage1, lineage2)
            
            if can_pair:
                if pair_type == "cross_major":
                    cross_major_combos.append((lineage1, lineage2))
                else:
                    within_major_combos.append((lineage1, lineage2))
    
    # Shuffle combinations for fairness
    random.shuffle(cross_major_combos)
    random.shuffle(within_major_combos)
    
    print(f"Valid cross-major combinations: {len(cross_major_combos)}")
    print(f"Valid within-major combinations: {len(within_major_combos)}")
    
    # Pass 1: Cross-major lineage pairs with balanced allocation
    print("\n=== PASS 1: Cross-major-lineage pairs (balanced) ===")
    cross_major_pairs = 0
    
    # Calculate how many pairs each combination should get
    if cross_major_combos:
        # Round-robin allocation to ensure fair distribution
        max_rounds = 50  # Prevent infinite loops
        round_num = 0
        
        while round_num < max_rounds:
            pairs_made_this_round = 0
            
            for lineage1, lineage2 in cross_major_combos:
                # Get unused samples from each lineage
                available1 = [s for s in lineage_groups[lineage1] if s['id'] not in used_samples]
                available2 = [s for s in lineage_groups[lineage2] if s['id'] not in used_samples]
                
                if available1 and available2:
                    # Make one pair from this combination
                    sample1 = random.choice(available1)
                    sample2 = random.choice(available2)
                    
                    pairs.append((sample1, sample2))
                    used_samples.add(sample1['id'])
                    used_samples.add(sample2['id'])
                    cross_major_pairs += 1
                    pairs_made_this_round += 1
                    
                    print(f"  Round {round_num + 1}: {lineage1} × {lineage2}")
            
            if pairs_made_this_round == 0:
                break  # No more pairs possible
            
            round_num += 1
    
    print(f"Cross-major-lineage pairs: {cross_major_pairs}")
    
    # Pass 2: Within-major lineage pairs
    print("\n=== PASS 2: Within-major-lineage pairs (balanced) ===")
    within_major_pairs = 0
    
    if within_major_combos:
        max_rounds = 50
        round_num = 0
        
        while round_num < max_rounds:
            pairs_made_this_round = 0
            
            for lineage1, lineage2 in within_major_combos:
                # Get unused samples from each lineage
                available1 = [s for s in lineage_groups[lineage1] if s['id'] not in used_samples]
                available2 = [s for s in lineage_groups[lineage2] if s['id'] not in used_samples]
                
                if available1 and available2:
                    # Make one pair from this combination
                    sample1 = random.choice(available1)
                    sample2 = random.choice(available2)
                    
                    pairs.append((sample1, sample2))
                    used_samples.add(sample1['id'])
                    used_samples.add(sample2['id'])
                    within_major_pairs += 1
                    pairs_made_this_round += 1
                    
                    print(f"  Round {round_num + 1}: {lineage1} × {lineage2}")
            
            if pairs_made_this_round == 0:
                break
            
            round_num += 1
    
    print(f"Within-major-lineage pairs: {within_major_pairs}")
    
    # Report final distribution
    print("\n=== FINAL PAIR DISTRIBUTION ===")
    pair_counts = {}
    for sample1, sample2 in pairs:
        combo = tuple(sorted([sample1['lineage'], sample2['lineage']]))
        pair_counts[combo] = pair_counts.get(combo, 0) + 1
    
    for combo, count in sorted(pair_counts.items()):
        print(f"  {combo[0]} × {combo[1]}: {count} pairs")
    
    # Report any forbidden pairs that were skipped
    print("\n=== FORBIDDEN PAIRS (same 2nd order sublineage) ===")
    forbidden_found = False
    for i in range(len(lineages)):
        for j in range(i + 1, len(lineages)):
            lineage1 = lineages[i]
            lineage2 = lineages[j]
            can_pair, pair_type = can_pair_lineages(lineage1, lineage2)
            if not can_pair:
                print(f"  SKIPPED: {lineage1} × {lineage2} ({pair_type})")
                forbidden_found = True
    
    if not forbidden_found:
        print("  None found")
    
    return pairs

def write_pairs(pairs, output_file):
    """Write pairs to output file"""
    with open(output_file, 'w') as f:
        f.write("Sample1_ID\tSample1_Mixed\tSample1_Lineage\tSample2_ID\tSample2_Mixed\tSample2_Lineage\tPair_Type\n")
        
        for sample1, sample2 in pairs:
            # Determine pair type
            can_pair, pair_type = can_pair_lineages(sample1['lineage'], sample2['lineage'])
            pair_type_label = "Cross-Major" if pair_type == "cross_major" else "Within-Major"
            
            f.write(f"{sample1['id']}\t{sample1['mixed_value']}\t{sample1['lineage']}\t")
            f.write(f"{sample2['id']}\t{sample2['mixed_value']}\t{sample2['lineage']}\t{pair_type_label}\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 pair_samples.py <input.tsv> [output.tsv] [max_pairs_per_lineage]")
        print("Example: python3 pair_samples.py samples.tsv pairs.tsv 5")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "sample_pairs.tsv"
    max_pairs = int(sys.argv[3]) if len(sys.argv) > 3 else None
    
    # Set random seed for reproducible results (optional)
    random.seed(42)
    
    print(f"Reading samples from {input_file}...")
    lineage_groups = read_samples(input_file)
    
    print(f"Creating pairs...")
    pairs = pair_samples(lineage_groups, max_pairs)
    
    print(f"\nCreated {len(pairs)} total pairs")
    print(f"Writing pairs to {output_file}...")
    write_pairs(pairs, output_file)
    
    print(f"Done! Check {output_file} for your sample pairs.")

if __name__ == "__main__":
    main()