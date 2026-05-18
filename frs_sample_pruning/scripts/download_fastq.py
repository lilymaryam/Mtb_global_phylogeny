#!/usr/bin/env python3

import sys
import subprocess
import requests

def get_run_accession(exp_id):
    """Convert experiment accession (SRX/ERX/DRX) to run accession using ENA API"""
    try:
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport"
        params = {
            'accession': exp_id,
            'result': 'read_run',
            'fields': 'run_accession'
        }
        
        response = requests.get(url, params=params, timeout=30)
        lines = response.text.strip().split('\n')
        
        # Skip header line, get first run accession
        if len(lines) >= 2:
            return lines[1].strip()
        
        return None
        
    except Exception as e:
        print(f"ERROR: Failed to lookup {exp_id}: {e}", file=sys.stderr)
        return None

def download_fastq(exp_id, output_dir, sample_name):
    """Download FASTQ files using fasterq-dump"""
    
    # Handle list format
    if exp_id.startswith('['):
        exp_id = exp_id.strip('[]\'\"').split(',')[0].strip('\'\"')
    
    run_id = get_run_accession(exp_id)
    if not run_id:
        print(f"ERROR: Could not convert {exp_id} to run accession", file=sys.stderr)
        return False
    
    cmd = [
        'fasterq-dump',
        '--split-files',
        '--force',
        '--outdir', output_dir,
        '--outfile', sample_name,
        run_id
    ]
    
    subprocess.run(cmd, check=True)
    return True

if __name__ == "__main__":
    exp_id = sys.argv[1]
    output_dir = sys.argv[2]
    sample_name = sys.argv[3]
    success = download_fastq(exp_id, output_dir, sample_name)
    sys.exit(0 if success else 1)