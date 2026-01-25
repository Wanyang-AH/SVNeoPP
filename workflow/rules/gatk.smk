"""WGS deduplication with GATK MarkDuplicatesSpark wrapper."""

gatk_index = wgs_index
GATK_IDS = WGS_IDS
GATK_PAIR_IDS = [gatk_index[sid]["pair_id"] for sid in GATK_IDS]

GATK_EXTRA = config["params"]["gatk"].get("extra", "")
GATK_JAVA_OPTS = config["params"]["gatk"].get("java_opts", "")
GATK_REF = config["references"]["gatk_fa"]
GATK_DICT = config["references"]["gatk_dict"]
GATK_SITES = config["references"]["gatk_sites"]


rule mark_duplicates_spark:
    input:
        bam="results/bwa/{pair_id}/{sample_id}.bam",
    output:
        bam="results/gatk/dedup/{pair_id}/{sample_id}.bam",
        metrics="results/gatk/dedup/{pair_id}/{sample_id}.metrics.txt",
    log:
        "logs/gatk/markduplicates/{pair_id}/{sample_id}.log",
    params:
        extra=GATK_EXTRA,
        java_opts=GATK_JAVA_OPTS,
    threads: 8
    resources:
        mem_mb=8000,
        tmpdir=TMPDIR,
    wrapper:
        "v7.6.0/bio/gatk/markduplicatesspark"


rule gatk_baserecalibrator:
    input:
        bam="results/gatk/dedup/{pair_id}/{sample_id}.bam",
        bai="results/gatk/dedup/{pair_id}/{sample_id}.bam.bai",
        ref=GATK_REF,
        dict=GATK_DICT,
        known=GATK_SITES,
    output:
        recal_table="results/gatk/bqsr/{pair_id}/{sample_id}.grp",
    log:
        "logs/gatk/baserecalibrator/{pair_id}/{sample_id}.log",
    params:
        extra=GATK_EXTRA,
        java_opts=GATK_JAVA_OPTS,
    threads: 4
    resources:
        mem_mb=4000,
        tmpdir=TMPDIR,
    wrapper:
        "v7.6.0/bio/gatk/baserecalibrator"


rule gatk_applybqsr:
    input:
        bam="results/gatk/dedup/{pair_id}/{sample_id}.bam",
        bai="results/gatk/dedup/{pair_id}/{sample_id}.bam.bai",
        ref=GATK_REF,
        dict=GATK_DICT,
        recal_table="results/gatk/bqsr/{pair_id}/{sample_id}.grp",
    output:
        bam="results/gatk/recal/{pair_id}/{sample_id}.bam",
    log:
        "logs/gatk/applybqsr/{pair_id}/{sample_id}.log",
    params:
        extra=GATK_EXTRA,
        java_opts=GATK_JAVA_OPTS,
        embed_ref=True,
    threads: 4
    resources:
        mem_mb=4000,
        tmpdir=TMPDIR,
    wrapper:
        "v9.0.0/bio/gatk/applybqsr"


rule gatk_all:
    input:
        expand(
            "results/gatk/dedup/{pair_id}/{sample_id}.bam",
            zip,
            pair_id=GATK_PAIR_IDS,
            sample_id=GATK_IDS,
        ),
        expand(
            "results/gatk/dedup/{pair_id}/{sample_id}.bam.bai",
            zip,
            pair_id=GATK_PAIR_IDS,
            sample_id=GATK_IDS,
        ),
        expand(
            "results/gatk/dedup/{pair_id}/{sample_id}.metrics.txt",
            zip,
            pair_id=GATK_PAIR_IDS,
            sample_id=GATK_IDS,
        ),
        expand(
            "results/gatk/bqsr/{pair_id}/{sample_id}.grp",
            zip,
            pair_id=GATK_PAIR_IDS,
            sample_id=GATK_IDS,
        ),
        expand(
            "results/gatk/recal/{pair_id}/{sample_id}.bam",
            zip,
            pair_id=GATK_PAIR_IDS,
            sample_id=GATK_IDS,
        ),
        expand(
            "results/gatk/recal/{pair_id}/{sample_id}.bam.bai",
            zip,
            pair_id=GATK_PAIR_IDS,
            sample_id=GATK_IDS,
        ),
