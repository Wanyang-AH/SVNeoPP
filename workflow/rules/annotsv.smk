"""SV annotation with AnnotSV."""

pairing = config.get("pairing", {})
pairs = list(pairing.keys())
ANNOTSV_BUILD = config["params"]["annotsv"].get("genome_build", "GRCh38")
ANNOTSV_EXTRA = config["params"]["annotsv"]["extra"]

rule annotsv:
    input:
        vcf="results/sv/{pair}/{pair}.svaba.sv.vcf",
    output:
        tsv="results/annotsv/{pair}/annotsv.tsv",
        vcf="results/annotsv/{pair}/annotsv.vcf",
    threads: 4
    resources:
        mem_mb=8000,
        tmpdir=TMPDIR,
    conda: "../envs/annotsv.yaml"
    shell:
        """
        mkdir -p results/annotsv/{wildcards.pair}
        AnnotSV \
            -SVinputFile {input.vcf} \
            -outputFile {output.tsv} \
            -genomeBuild {ANNOTSV_BUILD} \
            {ANNOTSV_EXTRA}
        cp {input.vcf} {output.vcf}
        """

rule annotsv_all:
    input:
        expand("results/annotsv/{pair}/annotsv.tsv", pair=pairs),
