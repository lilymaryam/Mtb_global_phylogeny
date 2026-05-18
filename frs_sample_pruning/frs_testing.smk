# Snakefile

configfile: "config.yaml"

DATA_DIR = config["data_dir"]
RAW_FASTQ_DIR = f"{DATA_DIR}/raw_fastq"
MIXED_FASTQ_DIR = f"{DATA_DIR}/mixed_fastq"
MYCO_RESULTS_DIR = f"{DATA_DIR}/myco_results"


# Load sample pairs
with open("sample_pairs_list.txt") as f:
    SAMPLE_PAIRS = [line.strip() for line in f if line.strip()]

# Define mixture ratios
#RATIOS = list(range(75, 100))
RATIOS = list(range(80, 100))
# TEST: Use only first 10 pairs and ratios 75-85
#SAMPLE_PAIRS = SAMPLE_PAIRS[:5]
#RATIOS = RATIOS[:10]

# Extract individual sample IDs
SAMPLES = set()
for pair in SAMPLE_PAIRS:
    sample1, sample2 = pair.split('_', 1)
    SAMPLES.add(sample1)
    SAMPLES.add(sample2)

SAMPLES = sorted(list(SAMPLES))

# Load mapping from TSV
SRX_MAPPING = {}
with open("biosample_srx_lookup.tsv") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            srx = parts[1].strip('[]\'\"').split(',')[0].strip('\'\"')
            SRX_MAPPING[parts[0]] = srx





print(f"Found {len(SAMPLE_PAIRS)} sample pairs")
print(f"Found {len(SAMPLES)} unique samples")
print(f"Ratios: {RATIOS[0]}-{RATIOS[-1]}")

'''
rule all:
    input:
        expand(f"{RAW_FASTQ_DIR}/{{sample}}_1.fastq", sample=SAMPLES),
        expand(f"{RAW_FASTQ_DIR}/{{sample}}_2.fastq", sample=SAMPLES),
        expand(f"{MIXED_FASTQ_DIR}/{{pair}}_{{ratio}}_R1.fastq", 
               pair=SAMPLE_PAIRS, ratio=RATIOS),
        expand(f"{MIXED_FASTQ_DIR}/{{pair}}_{{ratio}}_R2.fastq", 
               pair=SAMPLE_PAIRS, ratio=RATIOS),
        '''

rule all:
    input:
        expand(f"{MIXED_FASTQ_DIR}/{{pair}}_{{ratio}}_R1.fastq", 
               pair=SAMPLE_PAIRS, ratio=RATIOS),
        expand(f"{MIXED_FASTQ_DIR}/{{pair}}_{{ratio}}_R2.fastq", 
               pair=SAMPLE_PAIRS, ratio=RATIOS),
        expand(f"{MYCO_RESULTS_DIR}/{{pair}}_{{ratio}}/status.txt",
               pair=SAMPLE_PAIRS, ratio=RATIOS),
        f"{DATA_DIR}/frs_matrix.tsv",

rule download_fastq:
    output:
        r1=f"{RAW_FASTQ_DIR}/{{sample}}_1.fastq",
        r2=f"{RAW_FASTQ_DIR}/{{sample}}_2.fastq",
    log:
        "logs/download_{sample}.log"
    params:
        srx_id=lambda wildcards: SRX_MAPPING[wildcards.sample],
        output_dir=RAW_FASTQ_DIR,
    threads: 1
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=720,
        slurm_partition="medium",
        ncbi_api=1
    shell:
        """
        python3 scripts/download_fastq.py {params.srx_id} {params.output_dir} {wildcards.sample} > {log} 2>&1
        """

rule mix_fastq:
    input:
        r1_1=lambda wildcards: f"{RAW_FASTQ_DIR}/{wildcards.pair.split('_')[0]}_1.fastq",
        r2_1=lambda wildcards: f"{RAW_FASTQ_DIR}/{wildcards.pair.split('_')[0]}_2.fastq",
        r1_2=lambda wildcards: f"{RAW_FASTQ_DIR}/{wildcards.pair.split('_')[1]}_1.fastq",
        r2_2=lambda wildcards: f"{RAW_FASTQ_DIR}/{wildcards.pair.split('_')[1]}_2.fastq",
    output:
        r1=f"{MIXED_FASTQ_DIR}/{{pair}}_{{ratio}}_R1.fastq",
        r2=f"{MIXED_FASTQ_DIR}/{{pair}}_{{ratio}}_R2.fastq",
    resources:
        mem_mb=8000,
        cpus=8,
        runtime=720,
        slurm_partition="medium",
    log:
        "logs/mix_{pair}_{ratio}.log"
    shell:
        """
        python3 scripts/mix_fastqs.py {input.r1_1} {input.r2_1} {input.r1_2} {input.r2_2} \
            {output.r1} {output.r2} --ratio {wildcards.ratio} > {log} 2>&1
        """

rule make_myco_inputs:
    input:
        r1=f"{MIXED_FASTQ_DIR}/{{pair}}_{{ratio}}_R1.fastq",
        r2=f"{MIXED_FASTQ_DIR}/{{pair}}_{{ratio}}_R2.fastq",
    output:
        json=f"{MYCO_RESULTS_DIR}/{{pair}}_{{ratio}}/inputs.json",
    params:
        sample_name=lambda wildcards: f"{wildcards.pair}_{wildcards.ratio}",
    shell:
        """
        mkdir -p $(dirname {output.json})
        python3 scripts/make_myco_inputs.py {input.r1} {input.r2} {params.sample_name} {output.json}
        """

rule run_myco:
    input:
        json=f"{MYCO_RESULTS_DIR}/{{pair}}_{{ratio}}/inputs.json",
    output:
        marker=f"{MYCO_RESULTS_DIR}/{{pair}}_{{ratio}}/status.txt",
    log:
        "logs/myco_{pair}_{ratio}.log"
    params:
        output_dir=f"{MYCO_RESULTS_DIR}/{{pair}}_{{ratio}}",
        vcf=f"{MYCO_RESULTS_DIR}/{{pair}}_{{ratio}}/_LAST/out/tbd_vcfs/0/{{pair}}_{{ratio}}.vcf",
    resources:
        mem_mb=32000,
        cpus=16,
        runtime=720,
        slurm_partition="medium",
    shell:
        """
        miniwdl run myco_raw.wdl \
            --input {input.json} \
            --dir {params.output_dir} > {log} 2>&1 || true
        
        cd {params.output_dir}
        ls -dt *_myco 2>/dev/null | tail -n +2 | xargs rm -rf 2>/dev/null || true
        cd -
        
        if [[ -f {params.vcf} ]]; then
            echo "PASS" > {output.marker}
        else
            echo "FAILED" > {output.marker}
        fi
        """

rule calculate_frs:
    input:
        marker=f"{MYCO_RESULTS_DIR}/{{pair}}_{{ratio}}/status.txt",
    output:
        tsv=f"{DATA_DIR}/frs_results/{{pair}}_{{ratio}}.tsv",
    params:
        vcf=f"{MYCO_RESULTS_DIR}/{{pair}}_{{ratio}}/_LAST/out/tbd_vcfs/0/{{pair}}_{{ratio}}.vcf",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        
        status=$(cat {input.marker})
        
        if [[ "$status" == "PASS" ]] && [[ -f {params.vcf} ]]; then
            tot=$(grep -v "^#" {params.vcf} | wc -l)
            filt=$(bcftools filter -i 'FRS<.85' {params.vcf} 2>/dev/null | grep -v '^#' | wc -l)
            if [[ $tot -gt 0 ]]; then
                result=$(echo "scale=4; $filt/$tot" | bc)
            else
                result="NA"
            fi
            echo -e "{wildcards.pair}\t{wildcards.ratio}\t${{tot}}\t${{filt}}\t${{result}}" > {output.tsv}
        else
            echo -e "{wildcards.pair}\t{wildcards.ratio}\tNA\tNA\tNA" > {output.tsv}
        fi
        """

rule aggregate_frs:
    input:
        expand(f"{DATA_DIR}/frs_results/{{pair}}_{{ratio}}.tsv",
               pair=SAMPLE_PAIRS, ratio=RATIOS),
    output:
        matrix=f"{DATA_DIR}/frs_matrix.tsv",
    shell:
        """
        python3 scripts/build_frs_matrix.py {output.matrix} {input}
        """