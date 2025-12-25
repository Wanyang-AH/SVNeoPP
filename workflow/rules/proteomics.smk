"""Proteomics validation with FragPipe (placeholder stub)."""
import csv


def load_proteomics_samples():
    rows = []
    with open(config["samples"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("datatype") == "PROT":
                rows.append(row)
    return rows


prot_samples = load_proteomics_samples()
prot_ids = [row["sample_id"] for row in prot_samples]
prot_index = {row["sample_id"]: row for row in prot_samples}

rule fragpipe_search:
    input:
        raw=lambda wildcards: prot_index[wildcards.sample_id]["r1"],
        peptides=lambda wildcards: f"results/neoantigen/{prot_index[wildcards.sample_id]['pair_id']}/peptides.fasta",
    output:
        psm="results/proteomics/{sample_id}/psm.tsv",
    threads: 8
    resources:
        mem_mb=8000,
        tmpdir=TMPDIR,
    conda: "../envs/fragpipe.yaml"
    shell:
        """
        mkdir -p results/proteomics/{wildcards.sample_id}
        python - <<'PY'
import csv
from pathlib import Path
psm_path = Path("{output.psm}")
psm_path.parent.mkdir(parents=True, exist_ok=True)
with psm_path.open("w", newline="") as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(["peptide", "score", "source"])
    writer.writerow(["PLACEHOLDER", "0", "fragpipe_stub"])
PY
        """

rule merge_proteomics:
    input:
        psm="results/proteomics/{pair}_prot/psm.tsv",
        peptides="results/neoantigen/{pair}/peptides.tsv",
    output:
        merged="results/proteomics/{pair}/neoantigen_proteomics.tsv",
    threads: 2
    conda: "../envs/python.yaml"
    shell:
        """
        python scripts/merge_proteomics.py --psm {input.psm} --peptides {input.peptides} --output {output.merged}
        """

rule proteomics_all:
    input:
        expand("results/proteomics/{sample_id}/psm.tsv", sample_id=prot_ids),
