"""HLA typing with OptiType wrapper (tumor RNA only)."""

# 只使用 tumor RNA
optitype_index = tumor_rna_index
OPTITYPE_IDS = TUMOR_RNA_IDS
OPTITYPE_PAIR_IDS = [optitype_index[sid]["pair_id"] for sid in OPTITYPE_IDS]

OPTITYPE_EXTRA = config["params"]["hla"].get("optitype_extra", "")


rule optitype:
    input:
        reads=lambda wc: [
            f"results/fastp/{optitype_index[wc.sample_id]['pair_id']}/rna/{wc.sample_id}_R1.fastq.gz",
            f"results/fastp/{optitype_index[wc.sample_id]['pair_id']}/rna/{wc.sample_id}_R2.fastq.gz",
        ]
    output:
        pdf="results/hla/optitype/{pair_id}/{sample_id}_coverage_plot.pdf",
        tsv="results/hla/optitype/{pair_id}/{sample_id}_result.tsv",
    log:
        "logs/optitype/{pair_id}/{sample_id}.log"
    params:
        sequencing_type="rna",
        config="",
        extra=OPTITYPE_EXTRA,
    wrapper:
        "v2.9.1/bio/optitype"


rule hla_all:
    input:
        expand(
            "results/hla/optitype/{pair_id}/{sample_id}_result.tsv",
            zip,
            pair_id=OPTITYPE_PAIR_IDS,
            sample_id=OPTITYPE_IDS,
        ),
        expand(
            "results/hla/optitype/{pair_id}/{sample_id}_coverage_plot.pdf",
            zip,
            pair_id=OPTITYPE_PAIR_IDS,
            sample_id=OPTITYPE_IDS,
        )
