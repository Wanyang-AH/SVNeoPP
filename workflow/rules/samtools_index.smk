"""BAM indexing with samtools wrapper."""

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
