"""HLA typing with OptiType on RNA tumor samples."""
import csv


def load_hla_samples():
    rows = []
    with open(config["samples"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("datatype") == "RNA" and row.get("status") == "tumor":
                rows.append(row)
    return rows


hla_samples = load_hla_samples()
hla_ids = [row["sample_id"] for row in hla_samples]
hla_index = {row["sample_id"]: row for row in hla_samples}
OPTITYPE_EXTRA = config["params"]["hla"]["optitype_extra"]

rule optitype:
    input:
        r1="results/qc/{sample_id}_R1.fastq.gz",
        r2="results/qc/{sample_id}_R2.fastq.gz",
    output:
        result="results/hla/{sample_id}/OptiType_result.tsv",
    threads: 8
    resources:
        mem_mb=8000,
        tmpdir=TMPDIR,
    conda: "../envs/optitype.yaml"
    shell:
        """
        mkdir -p results/hla/{wildcards.sample_id}
        OptiTypePipeline.py \
            -i {input.r1} {input.r2} \
            -o results/hla/{wildcards.sample_id} \
            -p {wildcards.sample_id} \
            -v \
            {OPTITYPE_EXTRA}
        cp results/hla/{wildcards.sample_id}/{wildcards.sample_id}_result.tsv {output.result}
        """

rule hla_all:
    input:
        expand("results/hla/{sample_id}/OptiType_result.tsv", sample_id=hla_ids),
