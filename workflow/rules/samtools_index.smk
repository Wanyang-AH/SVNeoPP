"""BAM indexing with samtools wrapper."""

GATK_DEDUP_IDS = WGS_IDS
GATK_DEDUP_PAIR_IDS = [wgs_index[sid]["pair_id"] for sid in GATK_DEDUP_IDS]

rule samtools_index_wgs:
    input:
        "results/bwa/{pair_id}/{sample_id}.bam",
    output:
        "results/bwa/{pair_id}/{sample_id}.bam.bai",
    log:
        "logs/samtools_index/bwa/{pair_id}/{sample_id}.log",
    params:
        extra="",
    threads: 4
    wrapper:
        "v8.1.1/bio/samtools/index"


rule samtools_index_star:
    input:
        "results/star/{pair_id}/{sample_id}/aligned.bam",
    output:
        "results/star/{pair_id}/{sample_id}/aligned.bam.bai",
    log:
        "logs/samtools_index/star/{pair_id}/{sample_id}.log",
    params:
        extra="",
    threads: 4
    wrapper:
        "v8.1.1/bio/samtools/index"


rule samtools_index_gatk_dedup:
    input:
        "results/gatk/dedup/{pair_id}/{sample_id}.bam",
    output:
        "results/gatk/dedup/{pair_id}/{sample_id}.bam.bai",
    log:
        "logs/samtools_index/gatk/{pair_id}/{sample_id}.log",
    params:
        extra="",
    threads: 4
    wrapper:
        "v8.1.1/bio/samtools/index"


rule samtools_index_gatk_recal:
    input:
        "results/gatk/recal/{pair_id}/{sample_id}.bam",
    output:
        "results/gatk/recal/{pair_id}/{sample_id}.bam.bai",
    log:
        "logs/samtools_index/gatk_recal/{pair_id}/{sample_id}.log",
    params:
        extra="",
    threads: 4
    wrapper:
        "v8.1.1/bio/samtools/index"
