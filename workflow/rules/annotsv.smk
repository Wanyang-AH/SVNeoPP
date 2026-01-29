"""SV annotation with AnnotSV (SVaba output)."""

import os

ANNOTSV_SIF = "containers/annotsv.sif"
ANNOTSV_DB = config["references"]["annotsv_data"]

ANN_PAIR_IDS = sorted(set(SVABA_PAIR_IDS))
WORKDIR = os.path.abspath(".")


rule annotsv:
    input:
        vcf="results/svaba/{pair_id}/{pair_id}.svaba.somatic.sv.vcf.gz",
        db=ANNOTSV_DB,
    output:
        tsv="results/annotsv/{pair_id}/{pair_id}_AnnotSV.tsv",
    log:
        "logs/annotsv/{pair_id}.log",
    params:
        workdir=WORKDIR,
    threads: 1
    resources:
        mem_mb=2000,
        tmpdir=TMPDIR,
    shell:
        r"""
        mkdir -p results/annotsv/{wildcards.pair_id}
        singularity exec \
          -B {params.workdir}:/work \
          -B {input.db}:/database \
          --pwd /work/results/annotsv/{wildcards.pair_id} \
          {ANNOTSV_SIF} AnnotSV \
          -annotationsDir /database \
          -SVinputFile /work/{input.vcf} \
          -outputFile /work/{output.tsv} \
          > {log} 2>&1
        """


rule annotsv_all:
    input:
        expand(
            "results/annotsv/{pair_id}/{pair_id}_AnnotSV.tsv",
            pair_id=ANN_PAIR_IDS,
        )
