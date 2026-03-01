#!/usr/bin/env python3
"""Parse MHCflurry output and merge with peptide metadata."""

import argparse
import csv
import re
from pathlib import Path


def normalize_peptide(value: str | None) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", "", str(value).strip().upper())


def normalize_hla(value: str | None) -> str:
    text = (value or "").strip().replace(" ", "")
    if not text:
        return ""
    if text.upper().startswith("HLA-"):
        text = text[4:]
    if "*" not in text:
        m = re.match(r"^([A-Za-z]+)([0-9].*)$", text)
        if m:
            text = f"{m.group(1)}*{m.group(2)}"
    return f"HLA-{text.upper()}" if text else ""


def find_column(fieldnames: list[str], candidates: list[str]) -> str:
    lowered = {name.lower(): name for name in fieldnames}
    for cand in candidates:
        found = lowered.get(cand.lower())
        if found:
            return found
    raise ValueError(f"Missing required columns {candidates}; got fields: {fieldnames}")


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return fieldnames, rows


def build_map_index(path: Path) -> tuple[list[str], dict[tuple[str, str], dict[str, str]]]:
    fieldnames, rows = read_csv(path)
    allele_col = find_column(fieldnames, ["allele"])
    peptide_col = find_column(fieldnames, ["peptide", "peptide_sequence"])
    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (
            normalize_hla(row.get(allele_col)),
            normalize_peptide(row.get(peptide_col)),
        )
        if key[0] and key[1] and key not in index:
            index[key] = row
    return fieldnames, index


def run_pipeline(
    raw_csv: Path,
    map_csv: Path,
    out_csv: Path,
    pair_id: str,
    log_path: Path,
) -> None:
    raw_fields, raw_rows = read_csv(raw_csv)
    if not raw_fields:
        raise ValueError(f"MHCflurry raw output has no header: {raw_csv}")
    raw_allele_col = find_column(raw_fields, ["allele"])
    raw_peptide_col = find_column(raw_fields, ["peptide", "peptide_sequence"])

    map_fields, map_index = build_map_index(map_csv)
    skip_map_fields = {"pair_id", "allele", "peptide", "peptide_sequence"}
    extra_map_fields = [x for x in map_fields if x not in skip_map_fields and x not in raw_fields]
    out_fields = ["pair_id", "allele", "peptide_sequence"] + raw_fields + extra_map_fields

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    matched_rows = 0
    with out_csv.open("w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=out_fields)
        writer.writeheader()
        for raw_row in raw_rows:
            allele = normalize_hla(raw_row.get(raw_allele_col))
            peptide = normalize_peptide(raw_row.get(raw_peptide_col))
            map_row = map_index.get((allele, peptide), {})
            if map_row:
                matched_rows += 1

            out_row = {
                "pair_id": pair_id,
                "allele": allele,
                "peptide_sequence": peptide,
            }
            for field in raw_fields:
                out_row[field] = raw_row.get(field, "")
            for field in extra_map_fields:
                out_row[field] = map_row.get(field, "")
            writer.writerow(out_row)

    with log_path.open("w") as log:
        log.write("MHCflurry parse summary\n")
        log.write(f"pair_id: {pair_id}\n")
        log.write(f"raw_input: {raw_csv}\n")
        log.write(f"map_input: {map_csv}\n")
        log.write(f"raw_rows: {len(raw_rows)}\n")
        log.write(f"rows_with_metadata: {matched_rows}\n")
        log.write(f"allele_column: {raw_allele_col}\n")
        log.write(f"peptide_column: {raw_peptide_col}\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Parse MHCflurry predictions and merge with peptide metadata.")
    p.add_argument("--input-raw", required=True, help="Raw MHCflurry output CSV")
    p.add_argument("--input-map", required=True, help="Input map CSV from prepare step")
    p.add_argument("--output", required=True, help="Output merged CSV")
    p.add_argument("--pair-id", required=True, help="Pair ID")
    p.add_argument("--log", required=True, help="Log file")
    return p.parse_args()


def run_from_snakemake() -> None:
    pair_id = str(getattr(snakemake.wildcards, "pair_id", "unknown"))
    run_pipeline(
        raw_csv=Path(str(snakemake.input.raw_csv)),
        map_csv=Path(str(snakemake.input.map_csv)),
        out_csv=Path(str(snakemake.output.all_csv)),
        pair_id=pair_id,
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.all_csv)).with_suffix(".log"),
    )


def main() -> None:
    args = parse_args()
    run_pipeline(
        raw_csv=Path(args.input_raw),
        map_csv=Path(args.input_map),
        out_csv=Path(args.output),
        pair_id=args.pair_id,
        log_path=Path(args.log),
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
