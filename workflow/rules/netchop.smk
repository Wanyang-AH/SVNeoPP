"""Proteasomal cleavage prediction with NetChop."""

pairing = config.get("pairing", {})
pairs = list(pairing.keys())
NETCHOP_THRESHOLD = config["params"]["netchop"]["threshold"]
NETCHOP_EXTRA = config["params"]["netchop"]["extra"]

rule netchop:
    input:
        fasta="results/neoantigen/{pair}/peptides.fasta",
    output:
        tsv="results/netchop/{pair}/netchop.tsv",
    threads: 2
    resources:
        mem_mb=2000,
        tmpdir=TMPDIR,
    conda: "../envs/netchop.yaml"
    shell:
        """
        mkdir -p results/netchop/{wildcards.pair}
        netchop -i {input.fasta} -o {output.tsv} -t {NETCHOP_THRESHOLD} {NETCHOP_EXTRA}
        """

rule netchop_all:
    input:
        expand("results/netchop/{pair}/netchop.tsv", pair=pairs),
