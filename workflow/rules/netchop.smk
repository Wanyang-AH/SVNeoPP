"""Module A: convert annotsv2pep peptide CSV to NetChop FASTA input."""

import os

NETCHOP_PAIR_IDS = sorted(set(ANNOTSV2PEP_PAIR_IDS))
NETCHOP_HOME = os.path.abspath(config["references"]["netchop_home"])
NETCHOP_BIN = os.path.abspath(config["references"]["netchop_bin"])
NETCHOP_EXTRA = config["params"]["netchop"].get("extra", "")
NETCHOP_THRESHOLD = float(config["params"]["netchop"].get("threshold", 0.5))
NETCHOP_N_TERM_MIN = float(config["params"]["netchop"].get("n_term_min", NETCHOP_THRESHOLD))
NETCHOP_C_TERM_MIN = float(config["params"]["netchop"].get("c_term_min", NETCHOP_THRESHOLD))
NETCHOP_INTERNAL_MAX = float(config["params"]["netchop"].get("internal_max", NETCHOP_THRESHOLD))
NETCHOP_INTERNAL_COUNT_MAX = int(config["params"]["netchop"].get("internal_count_max", 0))


rule netchop_prepare_input:
    input:
        csv="results/annotsv2pep/{pair_id}/{pair_id}_AnnotSV.peptides.csv",
    output:
        fasta="results/netchop/{pair_id}/{pair_id}.input.fasta",
    log:
        "logs/netchop/{pair_id}.prepare.log",
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "scripts/netchop/netchop_prepare_input.py"


rule netchop_prepare_all:
    input:
        expand(
            "results/netchop/{pair_id}/{pair_id}.input.fasta",
            pair_id=NETCHOP_PAIR_IDS,
        )


rule netchop_prepare_l041:
    input:
        "results/netchop/l041/l041.input.fasta"


rule netchop_run:
    input:
        fasta="results/netchop/{pair_id}/{pair_id}.input.fasta",
    output:
        raw="results/netchop/{pair_id}/{pair_id}.raw.txt",
    log:
        "logs/netchop/{pair_id}.run.log",
    params:
        netchop_home=NETCHOP_HOME,
        netchop_bin=NETCHOP_BIN,
        extra=NETCHOP_EXTRA,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    shell:
        r"""
        if [ ! -x "{params.netchop_bin}" ]; then
          echo "ERROR: netChop binary not executable: {params.netchop_bin}" > {log}
          exit 1
        fi
        if [ ! -d "{params.netchop_home}" ]; then
          echo "ERROR: netChop home not found: {params.netchop_home}" > {log}
          exit 1
        fi
        "{params.netchop_bin}" -d "{params.netchop_home}" {params.extra} "{input.fasta}" > "{output.raw}" 2>> "{log}"
        """


rule netchop_run_all:
    input:
        expand(
            "results/netchop/{pair_id}/{pair_id}.raw.txt",
            pair_id=NETCHOP_PAIR_IDS,
        )


rule netchop_run_l041:
    input:
        "results/netchop/l041/l041.raw.txt"


rule netchop_parse:
    input:
        raw="results/netchop/{pair_id}/{pair_id}.raw.txt",
        fasta="results/netchop/{pair_id}/{pair_id}.input.fasta",
    output:
        metrics="results/netchop/{pair_id}/{pair_id}.metrics.csv",
    log:
        "logs/netchop/{pair_id}.parse.log",
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "scripts/netchop/netchop_parse.py"


rule netchop_filter:
    input:
        metrics="results/netchop/{pair_id}/{pair_id}.metrics.csv",
    output:
        filtered="results/netchop/{pair_id}/{pair_id}.filtered.csv",
    log:
        "logs/netchop/{pair_id}.filter.log",
    params:
        n_term_min=NETCHOP_N_TERM_MIN,
        c_term_min=NETCHOP_C_TERM_MIN,
        internal_max=NETCHOP_INTERNAL_MAX,
        internal_count_max=NETCHOP_INTERNAL_COUNT_MAX,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "scripts/netchop/netchop_filter.py"


rule netchop_filter_all:
    input:
        expand(
            "results/netchop/{pair_id}/{pair_id}.filtered.csv",
            pair_id=NETCHOP_PAIR_IDS,
        )


rule netchop_filter_l041:
    input:
        "results/netchop/l041/l041.filtered.csv"
