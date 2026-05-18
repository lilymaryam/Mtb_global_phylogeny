#!/usr/bin/env python3

import argparse
import subprocess
import random
import os
import sys


def count_reads(fastq_file):
    """Count reads in a FASTQ file"""
    with open(fastq_file) as f:
        return sum(1 for _ in f) // 4


def subsample_fastq(input_file, output_file, num_reads, seed):
    """Subsample a FASTQ file using seqtk"""
    cmd = ["seqtk", "sample", f"-s{seed}", input_file, str(num_reads)]
    with open(output_file, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


def concatenate_files(input_files, output_file):
    """Concatenate files into one"""
    with open(output_file, "w") as out:
        for f in input_files:
            with open(f) as inp:
                out.write(inp.read())


def mix_fastq_files(s1_r1, s1_r2, s2_r1, s2_r2, out_r1, out_r2, ratio, target_reads):
    """Mix two paired FASTQ samples at a specified ratio"""
    reads1 = min(count_reads(s1_r1), count_reads(s1_r2))
    reads2 = min(count_reads(s2_r1), count_reads(s2_r2))
    
    reads1_target = int((ratio / 100.0) * target_reads)
    reads2_target = int(((100 - ratio) / 100.0) * target_reads)
    
    # Adjust if we don't have enough reads
    if reads1_target > reads1 or reads2_target > reads2:
        max_possible = min(reads1, reads2) * 2
        reads1_target = int((ratio / 100.0) * max_possible)
        reads2_target = int(((100 - ratio) / 100.0) * max_possible)
    
    # Make sure output directory exists
    os.makedirs(os.path.dirname(out_r1), exist_ok=True)
    
    # Create temp files for subsampling
    temp_dir = os.path.dirname(out_r1)
    temp_s1_r1 = os.path.join(temp_dir, f"temp_s1_r1_{os.getpid()}.fastq")
    temp_s1_r2 = os.path.join(temp_dir, f"temp_s1_r2_{os.getpid()}.fastq")
    temp_s2_r1 = os.path.join(temp_dir, f"temp_s2_r1_{os.getpid()}.fastq")
    temp_s2_r2 = os.path.join(temp_dir, f"temp_s2_r2_{os.getpid()}.fastq")
    
    # Subsample (use same seed for paired files to keep pairs together)
    seed = random.randint(1, 10000)
    subsample_fastq(s1_r1, temp_s1_r1, reads1_target, seed)
    subsample_fastq(s1_r2, temp_s1_r2, reads1_target, seed)
    
    seed = random.randint(1, 10000)
    subsample_fastq(s2_r1, temp_s2_r1, reads2_target, seed)
    subsample_fastq(s2_r2, temp_s2_r2, reads2_target, seed)
    
    # Concatenate
    concatenate_files([temp_s1_r1, temp_s2_r1], out_r1)
    concatenate_files([temp_s1_r2, temp_s2_r2], out_r2)
    
    # Cleanup
    for f in [temp_s1_r1, temp_s1_r2, temp_s2_r1, temp_s2_r2]:
        os.remove(f)


def main():
    parser = argparse.ArgumentParser(description="Mix two paired FASTQ samples at a specified ratio")
    parser.add_argument("r1_1", help="Sample 1 R1")
    parser.add_argument("r2_1", help="Sample 1 R2")
    parser.add_argument("r1_2", help="Sample 2 R1")
    parser.add_argument("r2_2", help="Sample 2 R2")
    parser.add_argument("out_r1", help="Output R1")
    parser.add_argument("out_r2", help="Output R2")
    parser.add_argument("--ratio", type=int, required=True, help="Sample 1 percentage (0-100)")
    parser.add_argument("--target-reads", type=int, default=1_000_000, help="Target total reads")
    
    args = parser.parse_args()
    
    mix_fastq_files(
        args.r1_1, args.r2_1, args.r1_2, args.r2_2,
        args.out_r1, args.out_r2,
        args.ratio, args.target_reads,
    )


if __name__ == "__main__":
    main()