"""Immunogenicity scoring with DeepImmuno."""

pairing = config.get("pairing", {})
pairs = list(pairing.keys())

rule deepimmuno_score:
    input:
        binding="results/mhc/{pair}/binding_merged.tsv",
    output:
        tsv="results/immunogenicity/{pair}/deepimmuno.tsv",
    threads: 2
    conda: "../envs/deepimmuno.yaml"
    shell:
        """
        mkdir -p results/immunogenicity/{wildcards.pair}
        python - <<'PY'
import csv
from pathlib import Path
binding_path = Path("{input.binding}")
rows = []
if binding_path.exists():
    with binding_path.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            rows.append(row)

out_path = Path("{output.tsv}")
out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w", newline="") as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(["peptide", "allele", "score", "tool"])
    for row in rows:
        writer.writerow([row.get('peptide'), row.get('allele'), 0.5, 'deepimmuno_stub'])
PY
        """

rule deepimmuno_all:
    input:
        expand("results/immunogenicity/{pair}/deepimmuno.tsv", pair=pairs),
