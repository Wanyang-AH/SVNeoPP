"""Generate peptide candidates from annotated SVs and RNA expression."""

pairing = config.get("pairing", {})
pairs = list(pairing.keys())
PEPLENS = ",".join(str(x) for x in config["params"]["neoantigen"]["peptide_lengths"])
MIN_TPM = config["params"]["neoantigen"]["min_tpm"]

rule sv_to_peptides:
    input:
        annotsv="results/annotsv/{pair}/annotsv.tsv",
    output:
        tsv="results/neoantigen/{pair}/peptides.tsv",
        fasta="results/neoantigen/{pair}/peptides.fasta",
    threads: 2
    resources:
        mem_mb=2000,
        tmpdir=TMPDIR,
    conda: "../envs/python.yaml"
    shell:
        """
        python scripts/sv_to_peptides.py \
            --sv {input.annotsv} \
            --output-tsv {output.tsv} \
            --output-fasta {output.fasta} \
            --peptide-lengths {PEPLENS} \
            --pair-id {wildcards.pair} \
            --min-tpm {MIN_TPM}
        """

rule peptides_all:
    input:
        expand("results/neoantigen/{pair}/peptides.tsv", pair=pairs),
        expand("results/neoantigen/{pair}/peptides.fasta", pair=pairs),
