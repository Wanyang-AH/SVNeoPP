"""WGS alignment and BAM processing: BWA-MEM, sort, mark duplicates, BQSR."""
import csv


def load_wgs_samples():
    rows = []
    with open(config["samples"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("datatype") == "WGS":
                rows.append(row)
    return rows


wgs_samples = load_wgs_samples()
wgs_index = {row["sample_id"]: row for row in wgs_samples}
wgs_ids = [row["sample_id"] for row in wgs_samples]
BWA_THREADS = int(config["params"]["bwa"].get("threads", 8))
STAR_THREADS = int(config.get("threads", 4))
GENOME_FA = config["references"]["genome_fa"]
KNOWN_SITES = config["references"]["known_sites"]
SAMTOOLS_SORT_MEM = config["params"]["samtools"]["sort_memory"]
GATK_EXTRA = config["params"]["gatk"]["extra"]

rule bwa_mem_sort:
    input:
        r1="results/qc/{sample_id}_R1.fastq.gz",
        r2="results/qc/{sample_id}_R2.fastq.gz",
    output:
        bam=temp("results/wgs/{sample_id}.sorted.bam"),
    threads: BWA_THREADS
    resources:
        mem_mb=8000,
        tmpdir=TMPDIR,
    conda: "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} \
            -R '@RG\tID:{wildcards.sample_id}\tSM:{wildcards.sample_id}\tPL:ILLUMINA' \
            {GENOME_FA} {input.r1} {input.r2} | \
        samtools sort -@ {threads} -m {SAMTOOLS_SORT_MEM} -T {TMPDIR}/{wildcards.sample_id} -o {output.bam}
        """

rule mark_duplicates:
    input:
        bam="results/wgs/{sample_id}.sorted.bam",
    output:
        bam=temp("results/wgs/{sample_id}.dedup.bam"),
        bai=temp("results/wgs/{sample_id}.dedup.bam.bai"),
        metrics="results/wgs/{sample_id}.dedup.metrics.txt",
    threads: 6
    resources:
        mem_mb=8000,
        tmpdir=TMPDIR,
    conda: "../envs/gatk.yaml"
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --TMP_DIR {TMPDIR} \
            {GATK_EXTRA}
        samtools index -@ {threads} {output.bam}
        """

rule bqsr:
    input:
        bam="results/wgs/{sample_id}.dedup.bam",
        bai="results/wgs/{sample_id}.dedup.bam.bai",
    output:
        bam="results/wgs/{sample_id}.recal.bam",
        bai="results/wgs/{sample_id}.recal.bai",
        table=temp("results/wgs/{sample_id}.recal.table"),
    threads: 6
    resources:
        mem_mb=8000,
        tmpdir=TMPDIR,
    conda: "../envs/gatk.yaml"
    shell:
        """
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {GENOME_FA} \
            --known-sites {KNOWN_SITES} \
            -O {output.table} \
            {GATK_EXTRA} \
            --tmp-dir {TMPDIR}

        gatk ApplyBQSR \
            -I {input.bam} \
            -R {GENOME_FA} \
            --bqsr-recal-file {output.table} \
            -O {output.bam} \
            {GATK_EXTRA} \
            --tmp-dir {TMPDIR}
        samtools index -@ {threads} {output.bam}
        """

rule wgs_recal_all:
    input:
        expand("results/wgs/{sample_id}.recal.bam", sample_id=wgs_ids),
        expand("results/wgs/{sample_id}.recal.bai", sample_id=wgs_ids),
