"""Structural variant calling with SvABA (tumor-normal pairs)."""

pairing = config.get("pairing", {})
pairs = list(pairing.keys())
GENOME_FA = config["references"]["genome_fa"]
SVABA_EXTRA = config["params"]["svaba"]["extra"]

rule svaba_call:
    input:
        tumor_bam=lambda wildcards: f"results/wgs/{pairing[wildcards.pair]['tumor']}.recal.bam",
        tumor_bai=lambda wildcards: f"results/wgs/{pairing[wildcards.pair]['tumor']}.recal.bai",
        normal_bam=lambda wildcards: f"results/wgs/{pairing[wildcards.pair]['normal']}.recal.bam",
        normal_bai=lambda wildcards: f"results/wgs/{pairing[wildcards.pair]['normal']}.recal.bai",
    output:
        vcf="results/sv/{pair}/{pair}.svaba.sv.vcf",
        log="results/sv/{pair}/{pair}.svaba.log",
    threads: 12
    resources:
        mem_mb=16000,
        tmpdir=TMPDIR,
    conda: "../envs/svaba.yaml"
    shell:
        """
        mkdir -p results/sv/{wildcards.pair}
        cd results/sv/{wildcards.pair}
        svaba run \
            -t {input.tumor_bam} \
            -n {input.normal_bam} \
            -G {GENOME_FA} \
            -k {TMPDIR}/{wildcards.pair} \
            -a {wildcards.pair} \
            -p {threads} \
            {SVABA_EXTRA} \
            > {wildcards.pair}.svaba.log 2>&1
        """

rule sv_all:
    input:
        expand("results/sv/{pair}/{pair}.svaba.sv.vcf", pair=pairs),
