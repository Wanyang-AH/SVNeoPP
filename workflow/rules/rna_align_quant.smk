"""RNA alignment with STAR and quantification with kallisto."""
import csv


def load_rna_samples():
    rows = []
    with open(config["samples"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("datatype") == "RNA":
                rows.append(row)
    return rows


rna_samples = load_rna_samples()
rna_index = {row["sample_id"]: row for row in rna_samples}
rna_ids = [row["sample_id"] for row in rna_samples]

STAR_THREADS = int(config["params"]["star"].get("threads", 8))
STAR_INDEX = config["references"]["star_index"]
KALLISTO_INDEX = config["references"]["kallisto_index"]
KALLISTO_BOOT = config["params"]["kallisto"]["bootstrap"]
KALLISTO_EXTRA = config["params"]["kallisto"]["extra"]

rule star_align:
    input:
        r1="results/qc/{sample_id}_R1.fastq.gz",
        r2="results/qc/{sample_id}_R2.fastq.gz",
    output:
        bam="results/rna/{sample_id}.Aligned.sortedByCoord.out.bam",
        counts="results/rna/{sample_id}.ReadsPerGene.out.tab",
        sj="results/rna/{sample_id}.SJ.out.tab",
        log="results/rna/{sample_id}.Log.final.out",
    threads: STAR_THREADS
    resources:
        mem_mb=12000,
        tmpdir=TMPDIR,
    conda: "../envs/star.yaml"
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {STAR_INDEX} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix results/rna/{wildcards.sample_id}. \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            {config['params']['star']['extra']}
        """

rule kallisto_quant:
    input:
        r1="results/qc/{sample_id}_R1.fastq.gz",
        r2="results/qc/{sample_id}_R2.fastq.gz",
    output:
        abundance="results/rna/kallisto/{sample_id}/abundance.tsv",
        run_info="results/rna/kallisto/{sample_id}/run_info.json",
    threads: 8
    resources:
        mem_mb=4000,
        tmpdir=TMPDIR,
    conda: "../envs/kallisto.yaml"
    shell:
        """
        kallisto quant \
            -i {KALLISTO_INDEX} \
            -o results/rna/kallisto/{wildcards.sample_id} \
            -b {KALLISTO_BOOT} \
            --threads {threads} \
            {KALLISTO_EXTRA} \
            {input.r1} {input.r2}
        """

rule rna_all:
    input:
        expand("results/rna/{sample_id}.Aligned.sortedByCoord.out.bam", sample_id=rna_ids),
        expand("results/rna/kallisto/{sample_id}/abundance.tsv", sample_id=rna_ids),
