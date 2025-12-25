configfile: "config/config.yaml"

TMPDIR = config.get("tmpdir", "tmp")
FINAL_OUTPUTS = config.get("final_outputs", [])

rule all:
    input: FINAL_OUTPUTS

include: "workflow/rules/qc.smk"
include: "workflow/rules/wgs_align.smk"
include: "workflow/rules/rna_align_quant.smk"
include: "workflow/rules/hla.smk"
include: "workflow/rules/sv.smk"
include: "workflow/rules/annotsv.smk"
include: "workflow/rules/neoantigen.smk"
include: "workflow/rules/netchop.smk"
include: "workflow/rules/mhc_binding.smk"
include: "workflow/rules/deepimmuno.smk"
include: "workflow/rules/proteomics.smk"
include: "workflow/rules/report.smk"
include: "workflow/rules/tests.smk"
