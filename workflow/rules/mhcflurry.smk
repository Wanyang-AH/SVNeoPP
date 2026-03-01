"""MHCflurry integration: prepare input, ensure model, run prediction, parse/filter outputs."""

import os

MHCFLURRY_PAIR_IDS = sorted(set(NETCHOP_PAIR_IDS))
MHCFLURRY_LENGTHS = config["params"]["neoantigen"].get("peptide_lengths", [8, 9, 10, 11])
MHCFLURRY_RANK_THRESHOLD = float(config["params"]["mhcflurry"].get("rank_threshold", 2.0))
MHCFLURRY_EXTRA = config["params"]["mhcflurry"].get("extra", "")
MHCFLURRY_DOCKER_IMAGE = str(config["params"]["mhcflurry"].get("docker_image", "d365f27d62d2"))
MHCFLURRY_MODEL_NAME = str(config["params"]["mhcflurry"].get("model_name", "models_class1_presentation"))
MHCFLURRY_MODEL_INSTALL_DIR = os.path.abspath(
    os.path.expanduser(str(config["params"]["mhcflurry"].get("model_install_dir", "~/.local/share/mhcflurry")))
)
MHCFLURRY_MODEL_CHECK_SUBDIR = str(
    config["params"]["mhcflurry"].get("model_check_subdir", "4/2.2.0/models_class1_presentation")
)
MHCFLURRY_MODEL_TARBALL = os.path.abspath(
    config["references"].get(
        "mhcflurry_model_tarball",
        "resources/mhcflurry/models_class1_presentation.20200611.tar.bz2",
    )
)
MHCFLURRY_MODEL_STAMP = "results/mhcflurry/.model_ready.stamp"

MHCFLURRY_PAIR_TO_TUMOR_RNA = {}
for sample_id, sample_meta in tumor_rna_index.items():
    pair_id = sample_meta["pair_id"]
    if pair_id in MHCFLURRY_PAIR_TO_TUMOR_RNA and MHCFLURRY_PAIR_TO_TUMOR_RNA[pair_id] != sample_id:
        raise ValueError(
            f"Multiple tumor RNA samples mapped to pair_id={pair_id}: "
            f"{MHCFLURRY_PAIR_TO_TUMOR_RNA[pair_id]} and {sample_id}"
        )
    MHCFLURRY_PAIR_TO_TUMOR_RNA[pair_id] = sample_id


def get_mhcflurry_optitype_result_path(wildcards):
    pair_id = wildcards.pair_id
    sample_id = MHCFLURRY_PAIR_TO_TUMOR_RNA.get(pair_id)
    if sample_id is None:
        raise ValueError(f"No tumor RNA sample found for pair_id={pair_id} in config/samples.tsv")
    return f"results/hla/optitype/{pair_id}/{sample_id}_result.tsv"


rule mhcflurry_prepare_input:
    input:
        netchop_filtered="results/netchop/{pair_id}/{pair_id}.filtered.csv",
        optitype_tsv=get_mhcflurry_optitype_result_path,
    output:
        input_csv="results/mhcflurry/{pair_id}/{pair_id}.input.csv",
        map_csv="results/mhcflurry/{pair_id}/{pair_id}.input.map.csv",
    log:
        "logs/mhcflurry/{pair_id}.prepare.log",
    params:
        lengths=MHCFLURRY_LENGTHS,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/mhcflurry/prepare_input.py"


rule mhcflurry_model_ready:
    input:
        tarball=MHCFLURRY_MODEL_TARBALL,
    output:
        stamp=MHCFLURRY_MODEL_STAMP,
    log:
        "logs/mhcflurry/model_ready.log",
    params:
        docker_image=MHCFLURRY_DOCKER_IMAGE,
        model_name=MHCFLURRY_MODEL_NAME,
        model_install_dir=MHCFLURRY_MODEL_INSTALL_DIR,
        model_check_subdir=MHCFLURRY_MODEL_CHECK_SUBDIR,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.stamp}")" "$(dirname "{log}")" "{params.model_install_dir}"
        : > "{log}"
        if [ -d "{params.model_install_dir}/{params.model_check_subdir}" ]; then
          echo "MHCflurry model already installed: {params.model_install_dir}/{params.model_check_subdir}" >> "{log}"
          date -Is > "{output.stamp}"
          exit 0
        fi
        if [ ! -f "{input.tarball}" ]; then
          echo "ERROR: model tarball not found: {input.tarball}" >> "{log}"
          exit 1
        fi
        TARBALL_DIR="$(dirname "{input.tarball}")"
        echo "Installing MHCflurry model from tarball directory: $TARBALL_DIR" >> "{log}"
        docker run --rm -i \
          -e HOME=/home/docker \
          -v "{params.model_install_dir}:/home/docker/.local/share/mhcflurry" \
          -v "$TARBALL_DIR:/source_tarballs" \
          "{params.docker_image}" \
          mhcflurry-downloads fetch "{params.model_name}" --already-downloaded-dir /source_tarballs \
          >> "{log}" 2>&1
        if [ ! -d "{params.model_install_dir}/{params.model_check_subdir}" ]; then
          echo "ERROR: model installation check failed: {params.model_install_dir}/{params.model_check_subdir}" >> "{log}"
          exit 1
        fi
        date -Is > "{output.stamp}"
        """


rule mhcflurry_run:
    input:
        input_csv="results/mhcflurry/{pair_id}/{pair_id}.input.csv",
        model_stamp=MHCFLURRY_MODEL_STAMP,
    output:
        raw_csv="results/mhcflurry/{pair_id}/{pair_id}.raw.csv",
    log:
        "logs/mhcflurry/{pair_id}.run.log",
    params:
        docker_image=MHCFLURRY_DOCKER_IMAGE,
        model_install_dir=MHCFLURRY_MODEL_INSTALL_DIR,
        extra=MHCFLURRY_EXTRA,
    threads: 1
    resources:
        mem_mb=2000,
        tmpdir=TMPDIR,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.raw_csv}")" "$(dirname "{log}")"
        docker run --rm -i \
          -e HOME=/home/docker \
          -v "{params.model_install_dir}:/home/docker/.local/share/mhcflurry" \
          -v "$(pwd)":/work \
          -w /work \
          "{params.docker_image}" \
          mhcflurry-predict "{input.input_csv}" --out "{output.raw_csv}" {params.extra} \
          > "{log}" 2>&1
        """


rule mhcflurry_parse:
    input:
        raw_csv="results/mhcflurry/{pair_id}/{pair_id}.raw.csv",
        map_csv="results/mhcflurry/{pair_id}/{pair_id}.input.map.csv",
    output:
        all_csv="results/mhcflurry/{pair_id}/{pair_id}.all.csv",
    log:
        "logs/mhcflurry/{pair_id}.parse.log",
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/mhcflurry/parse.py"


rule mhcflurry_filter:
    input:
        all_csv="results/mhcflurry/{pair_id}/{pair_id}.all.csv",
    output:
        filtered="results/mhcflurry/{pair_id}/{pair_id}.filtered.csv",
    log:
        "logs/mhcflurry/{pair_id}.filter.log",
    params:
        rank_threshold=MHCFLURRY_RANK_THRESHOLD,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/mhcflurry/filter.py"


rule mhcflurry_all:
    input:
        expand(
            "results/mhcflurry/{pair_id}/{pair_id}.filtered.csv",
            pair_id=MHCFLURRY_PAIR_IDS,
        )


rule mhcflurry_filter_l041:
    input:
        "results/mhcflurry/l041/l041.filtered.csv"
