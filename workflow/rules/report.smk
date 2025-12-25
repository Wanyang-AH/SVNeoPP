"""Aggregation and reporting."""

pairing = config.get("pairing", {})
pairs = list(pairing.keys())

rule aggregate_report:
    input:
        peptides=expand("results/neoantigen/{pair}/peptides.tsv", pair=pairs),
        binding=expand("results/mhc/{pair}/binding_merged.tsv", pair=pairs),
        immunogenicity=expand("results/immunogenicity/{pair}/deepimmuno.tsv", pair=pairs),
        proteomics=expand("results/proteomics/{pair}/neoantigen_proteomics.tsv", pair=pairs),
    output:
        tsv="results/report/neoantigen_summary.tsv",
        xlsx="results/report/neoantigen_report.xlsx",
    params:
        pairs_str=" ".join(pairs),
        peptides_str=lambda wildcards, input: " ".join(input.peptides),
        binding_str=lambda wildcards, input: " ".join(input.binding),
        immuno_str=lambda wildcards, input: " ".join(input.immunogenicity),
        proteomics_str=lambda wildcards, input: " ".join(input.proteomics),
    threads: 2
    conda: "../envs/python.yaml"
    shell:
        """
        python scripts/aggregate_report.py \
            --pairs {params.pairs_str} \
            --peptides {params.peptides_str} \
            --binding {params.binding_str} \
            --immunogenicity {params.immuno_str} \
            --proteomics {params.proteomics_str} \
            --output-tsv {output.tsv} \
            --output-xlsx {output.xlsx}
        """

rule report_all:
    input:
        "results/report/neoantigen_summary.tsv",
        "results/report/neoantigen_report.xlsx",
