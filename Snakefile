configfile: "config/config.yaml"

# samples = workflow.source_path("config/samples.tsv")

TMPDIR = config.get("tmpdir", "tmp")
FINAL_OUTPUTS = config.get("final_outputs", [])



include: "workflow/rules/common.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/bwa.smk"
include: "workflow/rules/star.smk"
include: "workflow/rules/samtools_index.smk"
include: "workflow/rules/kallisto.smk"
include: "workflow/rules/hla.smk"
include: "workflow/rules/gatk.smk"
include: "workflow/rules/svaba.smk"


rule all:
    input: rules.kallisto_all.input,
           rules.star_all.input,
           rules.bwa_all.input,
           rules.hla_all.input,
           rules.gatk_all.input,
           rules.svaba_all.input,
