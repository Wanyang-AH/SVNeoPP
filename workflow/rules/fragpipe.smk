"""FragPipe integration: build custom FASTA, prepare manifest/workflow, and run headless search."""

import os


def _as_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


FRAGPIPE_PAIR_IDS = sorted({row["pair_id"] for row in prot_samples if row.get("pair_id")})
FRAGPIPE_WORKFLOW_DIR = os.path.abspath(
    config["references"].get("fragpipe_workflow_dir", "resources/fragpipe/workflows")
)
FRAGPIPE_REFERENCE_FASTA = os.path.abspath(
    config["references"].get(
        "fragpipe_reference_fasta",
        "resources/fragpipe/db/reference/uniprot-reviewed-UP000005640.fas",
    )
)
FRAGPIPE_DOCKER_IMAGE = str(config["params"]["fragpipe"].get("docker_image", "sibpt/fragpipe:anwy"))
FRAGPIPE_PHILOSOPHER_BIN = str(
    config["params"]["fragpipe"].get(
        "philosopher_bin_in_container",
        "/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2",
    )
)
FRAGPIPE_BIN = str(
    config["params"]["fragpipe"].get(
        "fragpipe_bin_in_container",
        "/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/bin/fragpipe",
    )
)
FRAGPIPE_TOOLS_FOLDER = str(config["params"]["fragpipe"].get("tools_folder_in_container", ""))
FRAGPIPE_RAM_GB = int(config["params"]["fragpipe"].get("ram_gb", 12))
FRAGPIPE_THREADS = int(config["params"]["fragpipe"].get("threads", 4))
FRAGPIPE_CONTAM = _as_bool(config["params"]["fragpipe"].get("contam", True))
FRAGPIPE_EXTRA = str(config["params"]["fragpipe"].get("extra", ""))


def get_fragpipe_workflow_template(wildcards):
    return os.path.join(FRAGPIPE_WORKFLOW_DIR, f"{wildcards.pair_id}.workflow")


def get_db_path_in_container(wildcards):
    pair_id = wildcards.pair_id
    return f"/work/results/fragpipe/{pair_id}/db/{pair_id}.decoys-contam.fas"


rule fragpipe_prepare_manifest:
    input:
        samples_tsv=config["samples"],
    output:
        manifest="results/fragpipe/{pair_id}/{pair_id}.fp-manifest",
    log:
        "logs/fragpipe/{pair_id}.manifest.log",
    params:
        mount_prefix="/work",
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/fragpipe/prepare_manifest.py"


rule fragpipe_prepare_custom_fasta:
    input:
        neo_fasta="results/netchop/{pair_id}/{pair_id}.input.fasta",
        reference_fasta=FRAGPIPE_REFERENCE_FASTA,
    output:
        custom_fasta="results/fragpipe/{pair_id}/db/{pair_id}.custom.fas",
    log:
        "logs/fragpipe/{pair_id}.custom_fasta.log",
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/fragpipe/prepare_custom_fasta.py"


rule fragpipe_build_db:
    input:
        custom_fasta="results/fragpipe/{pair_id}/db/{pair_id}.custom.fas",
    output:
        decoy_fasta="results/fragpipe/{pair_id}/db/{pair_id}.decoys-contam.fas",
        stamp="results/fragpipe/{pair_id}/db/.db_ready.stamp",
    log:
        "logs/fragpipe/{pair_id}.build_db.log",
    params:
        docker_image=FRAGPIPE_DOCKER_IMAGE,
        philosopher_bin=FRAGPIPE_PHILOSOPHER_BIN,
        contam=FRAGPIPE_CONTAM,
    threads: 1
    resources:
        mem_mb=2000,
        tmpdir=TMPDIR,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.decoy_fasta}")" "$(dirname "{output.stamp}")" "$(dirname "{log}")"
        REPO_ROOT="$(pwd)"
        WORKSPACE="${{REPO_ROOT}}/tmp/fragpipe/{wildcards.pair_id}/philosopher_workspace"
        mkdir -p "$WORKSPACE"
        rm -f "$WORKSPACE"/*decoy* "$WORKSPACE"/*decoys* "$WORKSPACE"/*contam* || true

        CONTAM_ARG=""
        if [ "{params.contam}" = "True" ] || [ "{params.contam}" = "true" ]; then
          CONTAM_ARG="--contam"
        fi

        {{
          echo "Workspace: $WORKSPACE"
          echo "Input custom FASTA: {input.custom_fasta}"
          echo "Output decoy FASTA: {output.decoy_fasta}"

          docker run --rm \
            -v "$WORKSPACE:/workspace" \
            -v "$REPO_ROOT:/work" \
            -w /workspace \
            "{params.docker_image}" \
            "{params.philosopher_bin}" workspace --init

          docker run --rm \
            -v "$WORKSPACE:/workspace" \
            -v "$REPO_ROOT:/work" \
            -w /workspace \
            "{params.docker_image}" \
            "{params.philosopher_bin}" database --custom "/work/{input.custom_fasta}" $CONTAM_ARG

          GENERATED="$(ls -1t "$WORKSPACE"/*decoy* "$WORKSPACE"/*decoys* 2>/dev/null | head -n1 || true)"
          if [ -z "$GENERATED" ]; then
            echo "ERROR: failed to find generated decoy database in $WORKSPACE"
            exit 1
          fi
          cp "$GENERATED" "{output.decoy_fasta}"

          docker run --rm \
            -v "$WORKSPACE:/workspace" \
            -w /workspace \
            "{params.docker_image}" \
            "{params.philosopher_bin}" workspace --clean || true
        }} > "{log}" 2>&1

        date -Is > "{output.stamp}"
        """


rule fragpipe_prepare_workflow:
    input:
        template=get_fragpipe_workflow_template,
        decoy_fasta="results/fragpipe/{pair_id}/db/{pair_id}.decoys-contam.fas",
    output:
        workflow="results/fragpipe/{pair_id}/{pair_id}.workflow",
    log:
        "logs/fragpipe/{pair_id}.workflow.log",
    params:
        db_path_in_container=get_db_path_in_container,
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir=TMPDIR,
    script:
        "../../scripts/fragpipe/prepare_workflow.py"


rule fragpipe_run:
    input:
        manifest="results/fragpipe/{pair_id}/{pair_id}.fp-manifest",
        workflow="results/fragpipe/{pair_id}/{pair_id}.workflow",
        decoy_fasta="results/fragpipe/{pair_id}/db/{pair_id}.decoys-contam.fas",
    output:
        done="results/fragpipe/{pair_id}/run/.done.stamp",
    log:
        "logs/fragpipe/{pair_id}.run.log",
    params:
        docker_image=FRAGPIPE_DOCKER_IMAGE,
        fragpipe_bin=FRAGPIPE_BIN,
        tools_folder=FRAGPIPE_TOOLS_FOLDER,
        ram_gb=FRAGPIPE_RAM_GB,
        extra=FRAGPIPE_EXTRA,
    threads: FRAGPIPE_THREADS
    resources:
        mem_mb=12000,
        tmpdir=TMPDIR,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.done}")" "$(dirname "{log}")"
        REPO_ROOT="$(pwd)"
        TOOLS_ARG=""
        if [ -n "{params.tools_folder}" ]; then
          TOOLS_ARG="--config-tools-folder {params.tools_folder}"
        fi

        docker run --rm \
          -v "$REPO_ROOT:/work" \
          -w /work \
          "{params.docker_image}" \
          "{params.fragpipe_bin}" --headless \
          --workflow "/work/{input.workflow}" \
          --manifest "/work/{input.manifest}" \
          --workdir "/work/results/fragpipe/{wildcards.pair_id}/run" \
          $TOOLS_ARG \
          --ram {params.ram_gb} \
          --threads {threads} \
          {params.extra} \
          > "{log}" 2>&1

        date -Is > "{output.done}"
        """


rule fragpipe_all:
    input:
        expand(
            "results/fragpipe/{pair_id}/run/.done.stamp",
            pair_id=FRAGPIPE_PAIR_IDS,
        )


rule fragpipe_l041:
    input:
        "results/fragpipe/l041/run/.done.stamp"
