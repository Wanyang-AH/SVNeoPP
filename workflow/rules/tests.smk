"""Validation hooks (lightweight sanity checks)."""

rule check_config:
    input:
        config["samples"],
    output:
        touch("results/validation/config_checked.txt"),
    threads: 1
    conda: "../envs/python.yaml"
    shell:
        """
        python - <<'PY'
import csv
import sys
from pathlib import Path
path = Path("{input}")
if not path.exists():
    sys.exit("samples.tsv missing")
with path.open() as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    rows = list(reader)
    if not rows:
        sys.exit("samples.tsv empty")
print(f"Loaded {len(rows)} samples from {path}")
PY
        """

rule validation_all:
    input:
        "results/validation/config_checked.txt",
