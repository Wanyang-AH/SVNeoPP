"""QC and trimming rules using fastp wrapper."""

rule fastp_pe:
    input:
        sample=lambda wildcards: [
            fastp_index[wildcards.sample_id]["r1"],
            fastp_index[wildcards.sample_id]["r2"],
        ]
    output:
        trimmed=[
            "results/qc/{sample_id}_R1.fastq.gz",
            "results/qc/{sample_id}_R2.fastq.gz",
        ],
        unpaired1="results/qc/{sample_id}.u1.fastq.gz",
        unpaired2="results/qc/{sample_id}.u2.fastq.gz",
        merged="results/qc/{sample_id}.merged.fastq.gz",
        failed="results/qc/{sample_id}.failed.fastq.gz",
        html="results/qc/{sample_id}.fastp.html",
        json="results/qc/{sample_id}.fastp.json",
    log:
        "results/qc/{sample_id}.fastp.log"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge",
    threads: 4
    resources:
        mem_mb=2000,
        tmpdir=TMPDIR,
    wrapper:
        "v7.1.0/bio/fastp"

rule fastp_all:
    input:
        expand("results/qc/{sample_id}_R1.fastq.gz", sample_id=FASTP_IDS),
        expand("results/qc/{sample_id}_R2.fastq.gz", sample_id=FASTP_IDS),
        expand("results/qc/{sample_id}.fastp.html", sample_id=FASTP_IDS),
        expand("results/qc/{sample_id}.fastp.json", sample_id=FASTP_IDS),
