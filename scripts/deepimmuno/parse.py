#!/usr/bin/env python3
"""Parse DeepImmuno raw TSV and merge with metadata map."""

import argparse
import csv
import re
from collections import defaultdict, deque
from pathlib import Path


def normalize_peptide(value: str | None) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", "", str(value).strip().upper())


def normalize_hla_for_deepimmuno(value: str | None) -> str:
    text = (value or "").strip().upper().replace(" ", "")
    if not text:
        return ""
    if text.startswith("HLA-"):
        text = text[4:]
    text = text.replace("*", "").replace(":", "")
    if len(text) < 5:
        return ""
    locus = text[0]
    digits = text[1:]
    if locus not in {"A", "B", "C"}:
        return ""
    if not digits.isdigit():
        return ""
    return f"HLA-{locus}*{digits}"


def find_column(fieldnames: list[str], candidates: list[str]) -> str:
    lowered = {name.lower(): name for name in fieldnames}
    for cand in candidates:
        found = lowered.get(cand.lower())
        if found:
            return found
    raise ValueError(f"Missing required columns {candidates}; got fields: {fieldnames}")


def read_map(path: Path) -> tuple[list[str], dict[tuple[str, str], deque[dict[str, str]]], int]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if "peptide_sequence" not in fieldnames or "deepimmuno_hla" not in fieldnames:
            raise ValueError(f"Map file must contain peptide_sequence and deepimmuno_hla: {path}")

        grouped: dict[tuple[str, str], deque[dict[str, str]]] = defaultdict(deque)
        total = 0
        for row in reader:
            key = (
                normalize_peptide(row.get("peptide_sequence")),
                normalize_hla_for_deepimmuno(row.get("deepimmuno_hla")),
            )
            grouped[key].append(row)
            total += 1
    return fieldnames, grouped, total


def read_raw(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    return fieldnames, rows


def run_pipeline(
    raw_tsv: Path,
    map_csv: Path,
    out_csv: Path,
    pair_id: str,
    log_path: Path,
) -> None:
    map_fieldnames, grouped_map, map_total = read_map(map_csv)
    raw_fieldnames, raw_rows = read_raw(raw_tsv)
    if not raw_fieldnames:
        raise ValueError(f"DeepImmuno raw output has no header: {raw_tsv}")

    pep_col = find_column(raw_fieldnames, ["peptide"])
    hla_col = find_column(raw_fieldnames, ["HLA", "hla"])
    score_col = find_column(raw_fieldnames, ["immunogenicity"])

    out_fields = map_fieldnames + ["deepimmuno_score"]
    out_rows: list[dict[str, str]] = []
    matched = 0
    unmatched_raw = 0

    for raw_row in raw_rows:
        peptide = normalize_peptide(raw_row.get(pep_col))
        hla = normalize_hla_for_deepimmuno(raw_row.get(hla_col))
        key = (peptide, hla)
        if grouped_map.get(key):
            map_row = grouped_map[key].popleft()
            matched += 1
        else:
            map_row = {
                "pair_id": pair_id,
                "peptide_sequence": peptide,
                "allele": "",
                "deepimmuno_hla": hla,
                "peptide_length": str(len(peptide)),
            }
            unmatched_raw += 1

        out_row = {field: map_row.get(field, "") for field in map_fieldnames}
        out_row["deepimmuno_score"] = (raw_row.get(score_col) or "").strip()
        out_rows.append(out_row)

    unmatched_map = 0
    for queue in grouped_map.values():
        while queue:
            map_row = queue.popleft()
            out_row = {field: map_row.get(field, "") for field in map_fieldnames}
            out_row["deepimmuno_score"] = ""
            out_rows.append(out_row)
            unmatched_map += 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields)
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        log.write("DeepImmuno parse summary\n")
        log.write(f"pair_id: {pair_id}\n")
        log.write(f"raw_input: {raw_tsv}\n")
        log.write(f"map_input: {map_csv}\n")
        log.write(f"map_rows: {map_total}\n")
        log.write(f"raw_rows: {len(raw_rows)}\n")
        log.write(f"matched_rows: {matched}\n")
        log.write(f"unmatched_raw_rows: {unmatched_raw}\n")
        log.write(f"unmatched_map_rows: {unmatched_map}\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Parse DeepImmuno output and merge with metadata.")
    p.add_argument("--input-raw", required=True, help="DeepImmuno raw TSV")
    p.add_argument("--input-map", required=True, help="Input metadata map CSV")
    p.add_argument("--output", required=True, help="Output merged CSV")
    p.add_argument("--pair-id", required=True, help="Pair ID")
    p.add_argument("--log", required=True, help="Log file path")
    return p.parse_args()


def run_from_snakemake() -> None:
    pair_id = str(getattr(snakemake.wildcards, "pair_id", "unknown"))
    run_pipeline(
        raw_tsv=Path(str(snakemake.input.raw_tsv)),
        map_csv=Path(str(snakemake.input.map_csv)),
        out_csv=Path(str(snakemake.output.all_csv)),
        pair_id=pair_id,
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.all_csv)).with_suffix(".log"),
    )


def main() -> None:
    args = parse_args()
    run_pipeline(
        raw_tsv=Path(args.input_raw),
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
