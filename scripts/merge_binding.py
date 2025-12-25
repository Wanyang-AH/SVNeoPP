#!/usr/bin/env python3
"""Merge NetMHCpan and MHCflurry prediction tables."""
import argparse
import csv
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Merge binding predictions")
    p.add_argument("--netmhc", required=True, help="NetMHCpan TSV")
    p.add_argument("--mhcflurry", required=True, help="MHCflurry TSV")
    p.add_argument("--output", required=True, help="Merged TSV output")
    return p.parse_args()


def read_rows(path):
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
    netmhc_rows = read_rows(args.netmhc)
    mhc_rows = read_rows(args.mhcflurry)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["peptide", "allele", "ic50", "rank", "tool"]
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in netmhc_rows + mhc_rows:
            writer.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
