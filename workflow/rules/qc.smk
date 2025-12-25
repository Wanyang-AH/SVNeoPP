"""QC and trimming rules using fastp."""
import csv


def load_samples(datatype):
    samples = []
    with open(config["samples"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("datatype") == datatype:
                samples.append(row)
    return samples


def build_index(rows):
    return {row["sample_id"]: row for row in rows}

wgs_samples = load_samples("WGS")
rna_samples = load_samples("RNA")
fastp_samples = wgs_samples + rna_samples
fastp_index = build_index(fastp_samples)
FASTP_IDS = [row["sample_id"] for row in fastp_samples]

rule fastp:
    input:
        r1=lambda wildcards: fastp_index[wildcards.sample_id]["r1"],
        r2=lambda wildcards: fastp_index[wildcards.sample_id]["r2"],
    params:
        extra=config["params"]["fastp"]["extra"],
    output:
        r1="results/qc/{sample_id}_R1.fastq.gz",
        r2="results/qc/{sample_id}_R2.fastq.gz",
        html="results/qc/{sample_id}.fastp.html",
        json="results/qc/{sample_id}.fastp.json",
    threads: 4
    resources:
        mem_mb=2000,
        tmpdir=TMPDIR,
    conda: "../envs/fastp.yaml"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --html {output.html} --json {output.json} \
            --thread {threads} {params.extra}
        """

rule fastp_all:
    input:
        expand("results/qc/{sample_id}_R1.fastq.gz", sample_id=FASTP_IDS),
        expand("results/qc/{sample_id}_R2.fastq.gz", sample_id=FASTP_IDS),
        expand("results/qc/{sample_id}.fastp.html", sample_id=FASTP_IDS),
        expand("results/qc/{sample_id}.fastp.json", sample_id=FASTP_IDS),
