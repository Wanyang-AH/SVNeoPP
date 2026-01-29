"""Somatic SV calling with SvABA (tumor-normal pairs)."""

import os

SVABA_SIF = "containers/svaba.sif"
SVABA_THREADS = config["params"]["svaba"].get("threads", 8)

SVABA_REF = config["references"]["gatk_fa"]
SVABA_DBSNP = config["references"]["gatk_dbsnp"]


def _to_ref_path(path):
    # Map resources/... to /ref/... inside the container
    return path.replace("resources/", "/ref/")


SVABA_TUMOR_BY_PAIR = {row["pair_id"]: sid for sid, row in tumor_wgs_index.items()}
SVABA_NORMAL_BY_PAIR = {row["pair_id"]: sid for sid, row in normal_wgs_index.items()}
SVABA_PAIR_IDS = sorted(set(SVABA_TUMOR_BY_PAIR) & set(SVABA_NORMAL_BY_PAIR))
WORKDIR = os.path.abspath(".")


rule svaba_run:
    input:
        tumor_bam=lambda wc: f"results/gatk/recal/{wc.pair_id}/{SVABA_TUMOR_BY_PAIR[wc.pair_id]}.bam",
        tumor_bai=lambda wc: f"results/gatk/recal/{wc.pair_id}/{SVABA_TUMOR_BY_PAIR[wc.pair_id]}.bam.bai",
        normal_bam=lambda wc: f"results/gatk/recal/{wc.pair_id}/{SVABA_NORMAL_BY_PAIR[wc.pair_id]}.bam",
        normal_bai=lambda wc: f"results/gatk/recal/{wc.pair_id}/{SVABA_NORMAL_BY_PAIR[wc.pair_id]}.bam.bai",
        ref=SVABA_REF,
        dbsnp=SVABA_DBSNP,
    output:
        sv_vcf="results/svaba/{pair_id}/{pair_id}.svaba.somatic.sv.vcf.gz",
    log:
        "logs/svaba/{pair_id}.log",
    params:
        ref=lambda wc: _to_ref_path(SVABA_REF), # ref/bwa_index/hg38.fa
        dbsnp=lambda wc: _to_ref_path(SVABA_DBSNP),
        workdir=WORKDIR,
    threads: SVABA_THREADS
    resources:
        mem_mb=16000,
        tmpdir=TMPDIR,
    shell:
        r"""
        mkdir -p results/svaba/{wildcards.pair_id}
        singularity exec \
          -B {params.workdir}/:/work \
          -B {params.workdir}/resources:/ref \
          --pwd /work/results/svaba/{wildcards.pair_id} \
          {SVABA_SIF} svaba run \
          -t /work/{input.tumor_bam} \
          -n /work/{input.normal_bam} \
          -G {params.ref} \
          -a {wildcards.pair_id} \
          -p {threads} \
          -z \
          -D {params.dbsnp} \
          > {log} 2>&1
        """


rule svaba_all:
    input:
        expand(
            "results/svaba/{pair_id}/{pair_id}.svaba.somatic.sv.vcf.gz",
            pair_id=SVABA_PAIR_IDS,
        )
