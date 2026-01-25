"""QC and trimming rules using fastp wrapper."""

# 源自common.smk文件
fastp_index = {**wgs_index, **rna_index}
FASTP_IDS = sorted(fastp_index.keys())
FASTP_PAIR_IDS = [fastp_index[sid]["pair_id"] for sid in FASTP_IDS]
FASTP_DTYPES = [fastp_index[sid]["datatype"] for sid in FASTP_IDS]

rule fastp_pe:
    input:
        sample=lambda wc: [
            fastp_index[wc.sample_id]["r1"],
            fastp_index[wc.sample_id]["r2"],
        ]
    output:
        trimmed=[
            "results/fastp/{pair_id}/{datatype}/{sample_id}_R1.fastq.gz",
            "results/fastp/{pair_id}/{datatype}/{sample_id}_R2.fastq.gz",
        ],
        unpaired1="results/fastp/{pair_id}/{datatype}/{sample_id}.u1.fastq.gz",
        unpaired2="results/fastp/{pair_id}/{datatype}/{sample_id}.u2.fastq.gz",
        merged="results/fastp/{pair_id}/{datatype}/{sample_id}.merged.fastq.gz",
        failed="results/fastp/{pair_id}/{datatype}/{sample_id}.failed.fastq.gz",
        html="results/fastp/{pair_id}/{datatype}/{sample_id}.fastp.html",
        json="results/fastp/{pair_id}/{datatype}/{sample_id}.fastp.json",
    log:
        "results/fastp/{pair_id}/{datatype}/{sample_id}.fastp.log"
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
        expand("results/fastp/{pair_id}/{datatype}/{sample_id}_R1.fastq.gz",
               zip, pair_id=FASTP_PAIR_IDS, datatype=FASTP_DTYPES, sample_id=FASTP_IDS),
        expand("results/fastp/{pair_id}/{datatype}/{sample_id}_R2.fastq.gz",
               zip, pair_id=FASTP_PAIR_IDS, datatype=FASTP_DTYPES, sample_id=FASTP_IDS),
        expand("results/fastp/{pair_id}/{datatype}/{sample_id}.fastp.html",
               zip, pair_id=FASTP_PAIR_IDS, datatype=FASTP_DTYPES, sample_id=FASTP_IDS),
        expand("results/fastp/{pair_id}/{datatype}/{sample_id}.fastp.json",
               zip, pair_id=FASTP_PAIR_IDS, datatype=FASTP_DTYPES, sample_id=FASTP_IDS),
