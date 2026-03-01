"""netMHCpan integration: prepare input, run prediction, parse/filter outputs."""

import os

NETMHCPAN_PAIR_IDS = sorted(set(NETCHOP_PAIR_IDS))
NETMHCPAN_HOME = os.path.abspath(config["references"]["netmhcpan_home"])
NETMHCPAN_BIN = os.path.abspath(config["references"]["netmhcpan_bin"])
NETMHCPAN_LENGTHS = config["params"]["netmhc"].get("peptide_lengths", [8, 9, 10, 11])
NETMHCPAN_LENGTHS_STR = ",".join(str(x) for x in NETMHCPAN_LENGTHS)
NETMHCPAN_RANK_EL_THRESHOLD = float(config["params"]["netmhc"].get("rank_el_threshold", 2.0))
NETMHCPAN_RANK_BA_THRESHOLD = float(config["params"]["netmhc"].get("rank_ba_threshold", 2.0))
NETMHCPAN_TMP_TEMPLATE = config["params"]["netmhc"].get("tmp_template", "/tmp/netMHCpanXXXXXX")
NETMHCPAN_EXTRA = config["params"]["netmhc"].get("extra", "")

PAIR_TO_TUMOR_RNA = {}
for sample_id, sample_meta in tumor_rna_index.items():
    pair_id = sample_meta["pair_id"]
    if pair_id in PAIR_TO_TUMOR_RNA and PAIR_TO_TUMOR_RNA[pair_id] != sample_id:
        raise ValueError(
            f"Multiple tumor RNA samples mapped to pair_id={pair_id}: "
            f"{PAIR_TO_TUMOR_RNA[pair_id]} and {sample_id}"
        )
    PAIR_TO_TUMOR_RNA[pair_id] = sample_id


def get_optitype_result_path(wildcards):
    pair_id = wildcards.pair_id
    sample_id = PAIR_TO_TUMOR_RNA.get(pair_id)
    if sample_id is None:
        raise ValueError(f"No tumor RNA sample found for pair_id={pair_id} in config/samples.tsv")
    return f"results/hla/optitype/{pair_id}/{sample_id}_result.tsv"


rule netmhcpan_prepare_input:
    input:
        netchop_filtered="results/netchop/{pair_id}/{pair_id}.filtered.csv",
        optitype_tsv=get_optitype_result_path,
    output:
        peptides="results/netmhcpan/{pair_id}/{pair_id}.input.pep",
        map_csv="results/netmhcpan/{pair_id}/{pair_id}.input.map.csv",
        alleles="results/netmhcpan/{pair_id}/{pair_id}.alleles.txt",
    log:
        "logs/netmhcpan/{pair_id}.prepare.log",
    params:
        lengths=NETMHCPAN_LENGTHS,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/netmhcpan/netmhcpan_prepare_input.py"


rule netmhcpan_run:
    input:
        peptides="results/netmhcpan/{pair_id}/{pair_id}.input.pep",
        alleles="results/netmhcpan/{pair_id}/{pair_id}.alleles.txt",
    output:
        raw="results/netmhcpan/{pair_id}/{pair_id}.raw.txt",
        xls="results/netmhcpan/{pair_id}/{pair_id}.xls",
    log:
        "logs/netmhcpan/{pair_id}.run.log",
    params:
        netmhcpan_home=NETMHCPAN_HOME,
        netmhcpan_bin=NETMHCPAN_BIN,
        lengths=NETMHCPAN_LENGTHS_STR,
        tmp_template=NETMHCPAN_TMP_TEMPLATE,
        extra=NETMHCPAN_EXTRA,
    threads: 1
    resources:
        mem_mb=2000,
        tmpdir=TMPDIR,
    shell:
        r"""
        if [ ! -x "{params.netmhcpan_bin}" ]; then
          echo "ERROR: netMHCpan binary not executable: {params.netmhcpan_bin}" > "{log}"
          exit 1
        fi
        if [ ! -d "{params.netmhcpan_home}" ]; then
          echo "ERROR: netMHCpan home not found: {params.netmhcpan_home}" > "{log}"
          exit 1
        fi
        ALLELES="$(tr -d '\r\n' < "{input.alleles}")"
        if [ -z "$ALLELES" ]; then
          echo "ERROR: empty allele list: {input.alleles}" > "{log}"
          exit 1
        fi
        NETMHCpan="{params.netmhcpan_home}" \
        "{params.netmhcpan_bin}" \
          -p "{input.peptides}" \
          -a "$ALLELES" \
          -l "{params.lengths}" \
          -BA \
          -rdir "{params.netmhcpan_home}" \
          -tdir "{params.tmp_template}" \
          {params.extra} \
          -xls \
          -xlsfile "{output.xls}" \
          > "{output.raw}" 2>> "{log}"
        """


rule netmhcpan_parse:
    input:
        xls="results/netmhcpan/{pair_id}/{pair_id}.xls",
        map_csv="results/netmhcpan/{pair_id}/{pair_id}.input.map.csv",
    output:
        all_csv="results/netmhcpan/{pair_id}/{pair_id}.all.csv",
    log:
        "logs/netmhcpan/{pair_id}.parse.log",
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/netmhcpan/netmhcpan_parse.py"


rule netmhcpan_filter:
    input:
        all_csv="results/netmhcpan/{pair_id}/{pair_id}.all.csv",
    output:
        filtered="results/netmhcpan/{pair_id}/{pair_id}.filtered.csv",
    log:
        "logs/netmhcpan/{pair_id}.filter.log",
    params:
        rank_el_threshold=NETMHCPAN_RANK_EL_THRESHOLD,
        rank_ba_threshold=NETMHCPAN_RANK_BA_THRESHOLD,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/netmhcpan/netmhcpan_filter.py"


rule netmhcpan_all:
    input:
        expand(
            "results/netmhcpan/{pair_id}/{pair_id}.filtered.csv",
            pair_id=NETMHCPAN_PAIR_IDS,
        )


rule netmhcpan_filter_l041:
    input:
        "results/netmhcpan/l041/l041.filtered.csv"
