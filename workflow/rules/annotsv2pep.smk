"""Generate candidate peptides from AnnotSV output."""

ANNOTSV2PEP_FA = config["references"]["gatk_fa"]
ANNOTSV2PEP_GTF = config["references"]["annotsv2pep_gtf"]
ANNOTSV2PEP_ID_MAP = config["references"]["annotsv2pep_id_map"]
ANNOTSV2PEP_DB = config["references"]["annotsv2pep_db"]
ANNOTSV2PEP_PAIR_IDS = sorted(set(ANN_PAIR_IDS))


rule annotsv2pep:
    input:
        tsv="results/annotsv/{pair_id}/{pair_id}_AnnotSV.tsv",
        fasta=ANNOTSV2PEP_FA,
        gtf=ANNOTSV2PEP_GTF,
        id_map=ANNOTSV2PEP_ID_MAP,
        db=ANNOTSV2PEP_DB,
    output:
        csv="results/annotsv2pep/{pair_id}/{pair_id}_AnnotSV.peptides.csv",
    log:
        "logs/annotsv2pep/{pair_id}.log",
    threads: 1
    resources:
        mem_mb=8000,
        tmpdir=TMPDIR,
    conda:
        "workflow/envs/annotsv2pep.yaml"
    script:
        "../../scripts/annotsv2pep.py"


rule annotsv2pep_all:
    input:
        expand(
            "results/annotsv2pep/{pair_id}/{pair_id}_AnnotSV.peptides.csv",
            pair_id=ANNOTSV2PEP_PAIR_IDS,
        )
