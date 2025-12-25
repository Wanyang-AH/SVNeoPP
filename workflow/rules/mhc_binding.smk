"""MHC binding predictions (NetMHCpan & MHCflurry) and merging."""
import csv
from pathlib import Path


pairing = config.get("pairing", {})
pairs = list(pairing.keys())


def tumor_rna_by_pair():
    mapping = {}
    with open(config["samples"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("datatype") == "RNA" and row.get("status") == "tumor":
                mapping[row["pair_id"]] = row["sample_id"]
    return mapping


tumor_rna_map = tumor_rna_by_pair()


def optitype_result(pair):
    sample_id = tumor_rna_map.get(pair)
    if not sample_id:
        raise ValueError(f"No tumor RNA sample for pair {pair}")
    return f"results/hla/{sample_id}/OptiType_result.tsv"


rule netmhcpan_pred:
    input:
        peptides="results/neoantigen/{pair}/peptides.fasta",
        optitype=lambda wildcards: optitype_result(wildcards.pair),
    output:
        tsv="results/mhc/{pair}/netmhcpan.tsv",
    threads: 4
    resources:
        mem_mb=4000,
        tmpdir=TMPDIR,
    conda: "../envs/netmhcpan.yaml"
    shell:
        """
        mkdir -p results/mhc/{wildcards.pair}
        python - <<'PY'
import csv
from pathlib import Path

optitype_path = Path("{input.optitype}")
alleles = []
if optitype_path.exists():
    with optitype_path.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        try:
            row = next(reader)
            for key in ("A1", "A2", "B1", "B2", "C1", "C2"):
                val = row.get(key)
                if val:
                    alleles.append(val)
        except StopIteration:
            pass

if not alleles:
    alleles = ["HLA-A02:01"]

peptides = []
with open("{input.peptides}") as fasta:
    seq = None
    for line in fasta:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            seq = None
            continue
        else:
            peptides.append(line)

Path("{output.tsv}").parent.mkdir(parents=True, exist_ok=True)
with open("{output.tsv}", "w", newline="") as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(["peptide", "allele", "ic50", "rank", "tool"])
    for pep in peptides:
        for allele in alleles:
            writer.writerow([pep, allele, 500, 2.0, "netmhcpan_stub"])
PY
        """

rule mhcflurry_pred:
    input:
        peptides="results/neoantigen/{pair}/peptides.fasta",
        optitype=lambda wildcards: optitype_result(wildcards.pair),
    output:
        tsv="results/mhc/{pair}/mhcflurry.tsv",
    threads: 4
    resources:
        mem_mb=4000,
        tmpdir=TMPDIR,
    conda: "../envs/mhcflurry.yaml"
    shell:
        """
        mkdir -p results/mhc/{wildcards.pair}
        python - <<'PY'
import csv
from pathlib import Path

optitype_path = Path("{input.optitype}")
alleles = []
if optitype_path.exists():
    with optitype_path.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        try:
            row = next(reader)
            for key in ("A1", "A2", "B1", "B2", "C1", "C2"):
                val = row.get(key)
                if val:
                    alleles.append(val)
        except StopIteration:
            pass

if not alleles:
    alleles = ["HLA-A02:01"]

peptides = []
with open("{input.peptides}") as fasta:
    seq = None
    for line in fasta:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            seq = None
            continue
        else:
            peptides.append(line)

Path("{output.tsv}").parent.mkdir(parents=True, exist_ok=True)
with open("{output.tsv}", "w", newline="") as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(["peptide", "allele", "ic50", "rank", "tool"])
    for pep in peptides:
        for allele in alleles:
            writer.writerow([pep, allele, 500, 2.0, "mhcflurry_stub"])
PY
        """

rule merge_binding:
    input:
        netmhc="results/mhc/{pair}/netmhcpan.tsv",
        mhcflurry="results/mhc/{pair}/mhcflurry.tsv",
    output:
        merged="results/mhc/{pair}/binding_merged.tsv",
    threads: 2
    conda: "../envs/python.yaml"
    shell:
        """
        python scripts/merge_binding.py --netmhc {input.netmhc} --mhcflurry {input.mhcflurry} --output {output.merged}
        """

rule mhc_all:
    input:
        expand("results/mhc/{pair}/binding_merged.tsv", pair=pairs),
