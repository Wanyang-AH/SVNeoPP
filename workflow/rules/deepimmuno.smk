"""DeepImmuno integration: prepare input, run prediction, parse/filter outputs."""

import os

DEEPIMMUNO_PAIR_IDS = sorted(set(NETMHCPAN_PAIR_IDS))
DEEPIMMUNO_HOME = os.path.abspath(config["references"].get("deepimmuno_home", "resources/DeepImmuno"))
DEEPIMMUNO_PYTHON = str(config["params"]["deepimmuno"].get("python", "python"))
DEEPIMMUNO_EXTRA = str(config["params"]["deepimmuno"].get("extra", ""))
DEEPIMMUNO_SUPPORTED_LENGTHS = config["params"]["deepimmuno"].get("supported_lengths", [9, 10])
DEEPIMMUNO_THRESHOLD = float(config["params"]["deepimmuno"].get("threshold", 0.0))


rule deepimmuno_prepare_input:
    input:
        netmhcpan_filtered="results/netmhcpan/{pair_id}/{pair_id}.filtered.csv",
    output:
        input_csv="results/deepimmuno/{pair_id}/{pair_id}.input.csv",
        map_csv="results/deepimmuno/{pair_id}/{pair_id}.input.map.csv",
    log:
        "logs/deepimmuno/{pair_id}.prepare.log",
    params:
        lengths=DEEPIMMUNO_SUPPORTED_LENGTHS,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/deepimmuno/prepare_input.py"


rule deepimmuno_run:
    input:
        input_csv="results/deepimmuno/{pair_id}/{pair_id}.input.csv",
    output:
        raw_tsv="results/deepimmuno/{pair_id}/{pair_id}.raw.tsv",
    log:
        "logs/deepimmuno/{pair_id}.run.log",
    params:
        deepimmuno_home=DEEPIMMUNO_HOME,
        python_bin=DEEPIMMUNO_PYTHON,
        extra=DEEPIMMUNO_EXTRA,
    threads: 1
    resources:
        mem_mb=8000,
        tmpdir=TMPDIR,
    conda:
        "workflow/envs/deepimmuno.yaml"
    script:
        "../../scripts/deepimmuno/run.py"


rule deepimmuno_parse:
    input:
        raw_tsv="results/deepimmuno/{pair_id}/{pair_id}.raw.tsv",
        map_csv="results/deepimmuno/{pair_id}/{pair_id}.input.map.csv",
    output:
        all_csv="results/deepimmuno/{pair_id}/{pair_id}.all.csv",
    log:
        "logs/deepimmuno/{pair_id}.parse.log",
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/deepimmuno/parse.py"


rule deepimmuno_filter:
    input:
        all_csv="results/deepimmuno/{pair_id}/{pair_id}.all.csv",
    output:
        filtered="results/deepimmuno/{pair_id}/{pair_id}.filtered.csv",
    log:
        "logs/deepimmuno/{pair_id}.filter.log",
    params:
        threshold=DEEPIMMUNO_THRESHOLD,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/deepimmuno/filter.py"


rule deepimmuno_all:
    input:
        expand(
            "results/deepimmuno/{pair_id}/{pair_id}.filtered.csv",
            pair_id=DEEPIMMUNO_PAIR_IDS,
        )


rule deepimmuno_filter_l041:
    input:
        "results/deepimmuno/l041/l041.filtered.csv"
