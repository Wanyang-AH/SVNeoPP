"""RNA quantification with kallisto wrapper."""

# 来源于 common.smk
rna_kallisto_index = rna_index
RNA_KALLISTO_IDS = RNA_IDS
RNA_PAIR_IDS = [rna_kallisto_index[sid]["pair_id"] for sid in RNA_KALLISTO_IDS]

KALLISTO_FASTA = config["references"]["kallisto_index"]
KALLISTO_INDEX = "resources/kallisto_index/kallisto.idx"
KALLISTO_EXTRA = config["params"]["kallisto"].get("extra", "")
KALLISTO_THREADS = int(config["params"]["kallisto"].get("threads", 1))


rule kallisto_index:
    input:
        fasta=KALLISTO_FASTA,
    output:
        index=KALLISTO_INDEX,
    params:
        extra="",
    log:
        "logs/kallisto_index.log",
    threads: 1
    wrapper:
        "v4.7.6/bio/kallisto/index"


rule kallisto_quant:
    input:
        fastq=lambda wc: [
            f"results/fastp/{rna_kallisto_index[wc.sample_id]['pair_id']}/rna/{wc.sample_id}_R1.fastq.gz",
            f"results/fastp/{rna_kallisto_index[wc.sample_id]['pair_id']}/rna/{wc.sample_id}_R2.fastq.gz",
        ],
        index=KALLISTO_INDEX,
    output:
        directory("results/kallisto/{pair_id}/{sample_id}"),
    log:
        "logs/kallisto_quant/{pair_id}/{sample_id}.log",
    params:
        extra=KALLISTO_EXTRA,
    threads: KALLISTO_THREADS
    wrapper:
        "v4.7.6/bio/kallisto/quant"


rule kallisto_all:
    input:
        expand(
            "results/kallisto/{pair_id}/{sample_id}",
            zip,
            pair_id=RNA_PAIR_IDS,
            sample_id=RNA_KALLISTO_IDS,
        )
