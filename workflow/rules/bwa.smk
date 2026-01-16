"""WGS alignment with BWA-MEM wrapper."""

# 来源于 common.smk
wgs_bwa_mem = wgs_index
WGS_BWA_IDS = WGS_IDS
WGS_PAIR_IDS = [wgs_bwa_mem[sid]["pair_id"] for sid in WGS_BWA_IDS]

rule bwa_mem:
    input:
        reads=lambda wc: [
            f"results/fastp/{wgs_bwa_mem[wc.sample_id]['pair_id']}/wgs/{wc.sample_id}_R1.fastq.gz",
            f"results/fastp/{wgs_bwa_mem[wc.sample_id]['pair_id']}/wgs/{wc.sample_id}_R2.fastq.gz",
        ],
        idx=multiext("resources/bwa_index/hg38.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/bwa/{pair_id}/{sample_id}.bam",
    log:
        "logs/bwa_mem/{pair_id}/{sample_id}.log",
    params:
        extra=r"-R '@RG\tID:{sample_id}\tSM:{sample_id}'",
        sorting="none",
        sort_order="queryname",
        sort_extra="",
    threads: 8
    wrapper:
        "v8.1.1/bio/bwa/mem"


rule bwa_all:
    input:
        expand(
            "results/bwa/{pair_id}/{sample_id}.bam",
            zip,
            pair_id=WGS_PAIR_IDS,
            sample_id=WGS_BWA_IDS,
        ),
