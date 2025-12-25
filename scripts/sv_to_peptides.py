#!/usr/bin/env python3
"""
Stub script to convert annotated SVs into candidate peptides.
This placeholder reads an AnnotSV TSV and emits a minimal peptide table/FASTA
based on provided peptide lengths; real junction/transcript logic should be
filled in later.
"""
import argparse
import csv
from pathlib import Path
from typing import List


def parse_args():
    parser = argparse.ArgumentParser(description="Generate peptide candidates from SVs")
    parser.add_argument("--sv", required=True, help="AnnotSV TSV input")
    parser.add_argument("--output-tsv", required=True, help="Output peptide table TSV")
    parser.add_argument("--output-fasta", required=True, help="Output peptide FASTA")
    parser.add_argument("--peptide-lengths", default="8,9,10,11", help="Comma-separated peptide lengths")
    parser.add_argument("--pair-id", required=True, help="Pair identifier")
    parser.add_argument("--min-tpm", type=float, default=0.0, help="Minimum TPM filter (placeholder)")
    return parser.parse_args()


def load_lengths(spec: str) -> List[int]:
    return [int(x) for x in spec.split(",") if x.strip()]


def read_sv_rows(path: Path):
    rows = []
    if not path.exists():
        return rows
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def main():
    args = parse_args()
    sv_rows = read_sv_rows(Path(args.sv))
    lengths = load_lengths(args.peptide_lengths)

    peptides = []
    for idx, sv in enumerate(sv_rows):
        sv_id = sv.get("AnnotSV ID", f"sv{idx+1}")
        for length in lengths:
            seq = f"PEPTIDE{length}_{idx+1}"[:length]
            peptides.append({
                "pair_id": args.pair_id,
                "sv_id": sv_id,
                "peptide": seq.ljust(length, "A"),
                "length": length,
                "source": "sv_stub"
            })

    # Ensure output directory exists
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_fasta).parent.mkdir(parents=True, exist_ok=True)

    with open(args.output_tsv, "w", newline="") as tsv_out:
        fieldnames = ["pair_id", "sv_id", "peptide", "length", "source"]
        writer = csv.DictWriter(tsv_out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in peptides:
            writer.writerow(row)

    with open(args.output_fasta, "w") as fasta_out:
        for row in peptides:
            fasta_out.write(f">{row['pair_id']}|{row['sv_id']}|len{row['length']}\n")
            fasta_out.write(f"{row['peptide']}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
