#!/usr/bin/env python3
"""Aggregate neoantigen evidence into summary tables (stub)."""
import argparse
import pandas as pd
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate neoantigen results")
    p.add_argument("--pairs", nargs="+", required=True, help="Pair identifiers")
    p.add_argument("--peptides", nargs="+", default=[], help="Peptide TSVs")
    p.add_argument("--binding", nargs="+", default=[], help="Binding merged TSVs")
    p.add_argument("--immunogenicity", nargs="+", default=[], help="DeepImmuno TSVs")
    p.add_argument("--proteomics", nargs="+", default=[], help="Proteomics evidence TSVs")
    p.add_argument("--output-tsv", required=True, help="Output summary TSV")
    p.add_argument("--output-xlsx", required=True, help="Output Excel report")
    return p.parse_args()


def read_peptides(peptide_paths):
    records = []
    for path in peptide_paths:
        p = Path(path)
        if not p.exists():
            continue
        pair = p.parts[-2] if len(p.parts) > 2 else "unknown"
        df = pd.read_csv(p, sep="\t")
        for _, row in df.iterrows():
            records.append({
                "pair": pair,
                "peptide": row.get("peptide", ""),
                "source": row.get("source", "sv"),
            })
    return records


def main():
    args = parse_args()
    peptide_records = read_peptides(args.peptides)
    if not peptide_records:
        peptide_records = [{"pair": pair, "peptide": "PLACEHOLDER", "source": "stub"} for pair in args.pairs]

    summary = pd.DataFrame(peptide_records)
    out_tsv = Path(args.output_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_tsv, sep="\t", index=False)

    out_xlsx = Path(args.output_xlsx)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_xlsx) as writer:
        summary.to_excel(writer, sheet_name="neoantigens", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
