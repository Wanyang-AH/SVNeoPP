#!/usr/bin/env python3
"""Prepare netMHCpan inputs from NetChop-filtered peptides and OptiType alleles."""

import argparse
import csv
import re
from pathlib import Path

VALID_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
HLA_COLUMNS = ("A1", "A2", "B1", "B2", "C1", "C2")


def normalize_peptide(value: str | None) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", "", str(value).strip().upper())


def normalize_hla(value: str | None) -> str:
    text = (value or "").strip().replace(" ", "")
    if not text:
        return ""
    if text.lower() in {"na", "nan", "none"}:
        return ""
    if text.upper().startswith("HLA-"):
        text = text[4:]
    text = text.replace("*", "")
    return f"HLA-{text}" if text else ""


def parse_lengths(raw: str | list[int] | tuple[int, ...]) -> set[int]:
    if isinstance(raw, (list, tuple)):
        return {int(x) for x in raw}
    text = str(raw).strip()
    if not text:
        return {8, 9, 10, 11}
    return {int(x.strip()) for x in text.split(",") if x.strip()}


def read_optitype_alleles(path: Path) -> list[str]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        first_row = next(reader, None)
        if first_row is None:
            raise ValueError(f"OptiType file is empty: {path}")

    alleles: list[str] = []
    seen: set[str] = set()
    for key in HLA_COLUMNS:
        allele = normalize_hla(first_row.get(key))
        if allele and allele not in seen:
            seen.add(allele)
            alleles.append(allele)
    if not alleles:
        raise ValueError(f"No valid HLA alleles found in {path}")
    return alleles


def prepare_inputs(
    in_csv: Path,
    optitype_tsv: Path,
    out_pep: Path,
    out_map_csv: Path,
    out_alleles: Path,
    pair_id: str,
    lengths: set[int],
    log_path: Path,
) -> None:
    out_pep.parent.mkdir(parents=True, exist_ok=True)
    out_map_csv.parent.mkdir(parents=True, exist_ok=True)
    out_alleles.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    alleles = read_optitype_alleles(optitype_tsv)

    total_rows = 0
    invalid_rows = 0
    skipped_length_rows = 0
    duplicate_rows = 0

    dedup_rows: list[dict[str, str]] = []
    seen_peptides: set[str] = set()

    with in_csv.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if "peptide_sequence" not in fieldnames:
            raise ValueError(f"Missing required column peptide_sequence in {in_csv}")

        for row in reader:
            total_rows += 1
            peptide = normalize_peptide(row.get("peptide_sequence"))
            if not peptide or not VALID_AA_RE.match(peptide):
                invalid_rows += 1
                continue
            if len(peptide) not in lengths:
                skipped_length_rows += 1
                continue
            if peptide in seen_peptides:
                duplicate_rows += 1
                continue
            seen_peptides.add(peptide)

            row = dict(row)
            row["pair_id"] = pair_id
            row["peptide_sequence"] = peptide
            row["peptide_length"] = str(len(peptide))
            dedup_rows.append(row)

    out_fields = ["pair_id", "peptide_sequence", "peptide_length"] + [
        x for x in (dedup_rows[0].keys() if dedup_rows else []) if x not in {"pair_id", "peptide_sequence", "peptide_length"}
    ]

    with out_pep.open("w") as pep_handle:
        for row in dedup_rows:
            pep_handle.write(f"{row['peptide_sequence']}\n")

    with out_map_csv.open("w", newline="") as map_handle:
        writer = csv.DictWriter(map_handle, fieldnames=out_fields)
        writer.writeheader()
        for row in dedup_rows:
            writer.writerow({k: row.get(k, "") for k in out_fields})

    with out_alleles.open("w") as allele_handle:
        allele_handle.write(",".join(alleles) + "\n")

    with log_path.open("w") as log:
        log.write("netMHCpan input preparation summary\n")
        log.write(f"input_csv: {in_csv}\n")
        log.write(f"optitype_tsv: {optitype_tsv}\n")
        log.write(f"pair_id: {pair_id}\n")
        log.write(f"lengths: {sorted(lengths)}\n")
        log.write(f"total_rows: {total_rows}\n")
        log.write(f"written_peptides: {len(dedup_rows)}\n")
        log.write(f"invalid_rows: {invalid_rows}\n")
        log.write(f"skipped_length_rows: {skipped_length_rows}\n")
        log.write(f"duplicate_rows: {duplicate_rows}\n")
        log.write(f"alleles: {','.join(alleles)}\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare netMHCpan input files from NetChop and OptiType outputs.")
    p.add_argument("--input-csv", required=True, help="Input NetChop filtered CSV")
    p.add_argument("--optitype", required=True, help="OptiType result TSV")
    p.add_argument("--output-pep", required=True, help="Output peptide list (.pep)")
    p.add_argument("--output-map", required=True, help="Output peptide metadata map CSV")
    p.add_argument("--output-alleles", required=True, help="Output allele list file")
    p.add_argument("--pair-id", required=True, help="Pair ID")
    p.add_argument("--lengths", default="8,9,10,11", help="Comma-separated peptide lengths")
    p.add_argument("--log", required=True, help="Log file path")
    return p.parse_args()


def run_from_snakemake() -> None:
    pair_id = str(getattr(snakemake.wildcards, "pair_id", "unknown"))
    prepare_inputs(
        in_csv=Path(str(snakemake.input.netchop_filtered)),
        optitype_tsv=Path(str(snakemake.input.optitype_tsv)),
        out_pep=Path(str(snakemake.output.peptides)),
        out_map_csv=Path(str(snakemake.output.map_csv)),
        out_alleles=Path(str(snakemake.output.alleles)),
        pair_id=pair_id,
        lengths=parse_lengths(snakemake.params.lengths),
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.peptides)).with_suffix(".log"),
    )


def main() -> None:
    args = parse_args()
    prepare_inputs(
        in_csv=Path(args.input_csv),
        optitype_tsv=Path(args.optitype),
        out_pep=Path(args.output_pep),
        out_map_csv=Path(args.output_map),
        out_alleles=Path(args.output_alleles),
        pair_id=args.pair_id,
        lengths=parse_lengths(args.lengths),
        log_path=Path(args.log),
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
