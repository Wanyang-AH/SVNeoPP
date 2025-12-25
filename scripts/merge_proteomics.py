#!/usr/bin/env python3
"""Merge proteomics evidence with neoantigen peptides."""
import argparse
import csv
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Merge proteomics evidence")
    p.add_argument("--psm", required=True, help="Proteomics PSM/peptide table")
    p.add_argument("--peptides", required=True, help="Neoantigen peptide table")
    p.add_argument("--output", required=True, help="Merged TSV output")
    return p.parse_args()


def read_tsv(path):
    rows = []
    if not Path(path).exists():
        return rows
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def main():
    args = parse_args()
    psm_rows = read_tsv(args.psm)
    peptide_rows = read_tsv(args.peptides)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["peptide", "evidence", "source"]
    with out_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        peptide_set = {row.get("peptide") for row in peptide_rows}
        for pep in peptide_set:
            writer.writerow({"peptide": pep, "evidence": "placeholder", "source": "fragpipe_stub"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
