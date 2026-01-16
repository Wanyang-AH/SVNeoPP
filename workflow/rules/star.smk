"""RNA alignment with STAR wrapper."""

# 来源于 common.smk
rna_star_index = rna_index
RNA_STAR_IDS = RNA_IDS
RNA_PAIR_IDS = [rna_star_index[sid]["pair_id"] for sid in RNA_STAR_IDS]

STAR_INDEX = config["references"]["star_index"]
STAR_THREADS = int(config["params"]["star"].get("threads", 8))
STAR_EXTRA = config["params"]["star"].get("extra", "")


rule star_pe_multi:
    input:
        fq1=lambda wc: [
            f"results/fastp/{rna_star_index[wc.sample_id]['pair_id']}/rna/{wc.sample_id}_R1.fastq.gz"
        ],
        fq2=lambda wc: [
            f"results/fastp/{rna_star_index[wc.sample_id]['pair_id']}/rna/{wc.sample_id}_R2.fastq.gz"
        ],
        idx=STAR_INDEX,
    output:
        aln="results/star/{pair_id}/{sample_id}/aligned.sam",
        log="results/star/{pair_id}/{sample_id}/Log.out",
        sj="results/star/{pair_id}/{sample_id}/SJ.out.tab",
        unmapped=[
            "results/star/{pair_id}/{sample_id}/unmapped.1.fastq.gz",
            "results/star/{pair_id}/{sample_id}/unmapped.2.fastq.gz",
        ],
    log:
        "logs/star/{pair_id}/{sample_id}.log",
    params:
        extra=STAR_EXTRA,
    threads: STAR_THREADS
    wrapper:
        "v3.3.7/bio/star/align"


rule star_all:
    input:
        expand(
            "results/star/{pair_id}/{sample_id}/aligned.sam",
            zip,
            pair_id=RNA_PAIR_IDS,
            sample_id=RNA_STAR_IDS,
        ),
        expand(
            "results/star/{pair_id}/{sample_id}/SJ.out.tab",
            zip,
            pair_id=RNA_PAIR_IDS,
            sample_id=RNA_STAR_IDS,
        ),
