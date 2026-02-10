#!/usr/bin/env python3
"""Prepare NetChop FASTA input from annotsv2pep peptide CSV."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable

VALID_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def clean_value(value: str | None, default: str = "NA") -> str:
    if value is None:
        return default
    text = str(value).strip()
    if not text:
        return default
    return text.replace("|", "/")


def normalize_peptide(value: str | None) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", "", str(value).strip().upper())


def read_csv_rows(csv_path: Path) -> Iterable[dict[str, str]]:
    with csv_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            yield row


def write_fasta(in_csv: Path, out_fasta: Path, sample_name: str, log_path: Path) -> None:
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    invalid_rows = 0
    duplicate_rows = 0
    written_rows = 0
    seen = set()

    with out_fasta.open("w") as fasta_handle:
        for row in read_csv_rows(in_csv):
            total_rows += 1
            peptide = normalize_peptide(row.get("peptide"))
            if not peptide or not VALID_AA_RE.match(peptide):
                invalid_rows += 1
                continue
            if peptide in seen:
                duplicate_rows += 1
                continue
            seen.add(peptide)

            sv_id = clean_value(row.get("sv_id"))
            gene_name = clean_value(row.get("gene_name"))
            pep_source = clean_value(row.get("peptide_source"))
            sv_type = clean_value(row.get("sv_type"))
            idx = written_rows + 1
            header = f"{sv_id}|{gene_name}|{pep_source}|{sv_type}|{sample_name}|idx{idx}"

            fasta_handle.write(f">{header}\n")
            fasta_handle.write(f"{peptide}\n")
            written_rows += 1

    with log_path.open("w") as log_handle:
        log_handle.write("NetChop input preparation summary\n")
        log_handle.write(f"input_csv: {in_csv}\n")
        log_handle.write(f"output_fasta: {out_fasta}\n")
        log_handle.write(f"sample: {sample_name}\n")
        log_handle.write(f"total_rows: {total_rows}\n")
        log_handle.write(f"written_rows: {written_rows}\n")
        log_handle.write(f"duplicate_rows: {duplicate_rows}\n")
        log_handle.write(f"invalid_rows: {invalid_rows}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare FASTA for NetChop from peptide CSV.")
    parser.add_argument("--input-csv", required=True, help="Input peptide CSV")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA path")
    parser.add_argument("--sample", required=True, help="Sample or pair identifier")
    parser.add_argument("--log", required=True, help="Log file path")
    return parser.parse_args()


def run_from_snakemake() -> None:
    input_csv = Path(str(snakemake.input.csv))
    output_fasta = Path(str(snakemake.output.fasta))
    log_path = Path(str(snakemake.log[0])) if snakemake.log else output_fasta.with_suffix(".log")
    sample = str(getattr(snakemake.wildcards, "pair_id", output_fasta.stem))
    write_fasta(input_csv, output_fasta, sample, log_path)


def main() -> None:
    args = parse_args()
    write_fasta(
        in_csv=Path(args.input_csv),
        out_fasta=Path(args.output_fasta),
        sample_name=args.sample,
        log_path=Path(args.log),
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
