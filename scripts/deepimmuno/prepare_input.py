#!/usr/bin/env python3
"""Prepare DeepImmuno input from netMHCpan filtered output."""

import argparse
import csv
import re
from pathlib import Path

VALID_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYX]+$")


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


def parse_lengths(raw: str | list[int] | tuple[int, ...]) -> set[int]:
    if isinstance(raw, (list, tuple)):
        return {int(x) for x in raw}
    text = str(raw).strip()
    if not text:
        return {9, 10}
    return {int(x.strip()) for x in text.split(",") if x.strip()}


def prepare_inputs(
    in_csv: Path,
    out_input_csv: Path,
    out_map_csv: Path,
    pair_id: str,
    lengths: set[int],
    log_path: Path,
) -> None:
    out_input_csv.parent.mkdir(parents=True, exist_ok=True)
    out_map_csv.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    invalid_rows = 0
    skipped_length_rows = 0
    duplicate_rows = 0
    input_fieldnames: list[str] = []
    prepared_rows: list[dict[str, str]] = []
    seen_pairs: set[tuple[str, str]] = set()

    with in_csv.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        input_fieldnames = reader.fieldnames or []
        if "peptide_sequence" not in input_fieldnames:
            raise ValueError(f"Missing required column peptide_sequence in {in_csv}")
        if "allele" not in input_fieldnames:
            raise ValueError(f"Missing required column allele in {in_csv}")

        for row in reader:
            total_rows += 1
            peptide = normalize_peptide(row.get("peptide_sequence"))
            allele = normalize_hla_for_deepimmuno(row.get("allele"))
            if not peptide or not allele or not VALID_AA_RE.match(peptide):
                invalid_rows += 1
                continue
            if len(peptide) not in lengths:
                skipped_length_rows += 1
                continue
            key = (peptide, allele)
            if key in seen_pairs:
                duplicate_rows += 1
                continue
            seen_pairs.add(key)

            item = dict(row)
            item["pair_id"] = pair_id
            item["peptide_sequence"] = peptide
            item["peptide_length"] = str(len(peptide))
            item["deepimmuno_hla"] = allele
            prepared_rows.append(item)

    with out_input_csv.open("w", newline="") as out_csv:
        for row in prepared_rows:
            out_csv.write(f"{row['peptide_sequence']},{row['deepimmuno_hla']}\n")

    base_fields = ["pair_id", "peptide_sequence", "allele", "deepimmuno_hla", "peptide_length"]
    skip_fields = set(base_fields)
    extra_fields = [x for x in input_fieldnames if x not in skip_fields]
    out_fields = base_fields + extra_fields

    with out_map_csv.open("w", newline="") as out_map:
        writer = csv.DictWriter(out_map, fieldnames=out_fields)
        writer.writeheader()
        for row in prepared_rows:
            out_row = {
                "pair_id": pair_id,
                "peptide_sequence": row["peptide_sequence"],
                "allele": row.get("allele", ""),
                "deepimmuno_hla": row["deepimmuno_hla"],
                "peptide_length": row["peptide_length"],
            }
            for field in extra_fields:
                out_row[field] = row.get(field, "")
            writer.writerow(out_row)

    with log_path.open("w") as log:
        log.write("DeepImmuno input preparation summary\n")
        log.write(f"input_csv: {in_csv}\n")
        log.write(f"pair_id: {pair_id}\n")
        log.write(f"lengths: {sorted(lengths)}\n")
        log.write(f"total_rows: {total_rows}\n")
        log.write(f"written_rows: {len(prepared_rows)}\n")
        log.write(f"invalid_rows: {invalid_rows}\n")
        log.write(f"skipped_length_rows: {skipped_length_rows}\n")
        log.write(f"duplicate_rows: {duplicate_rows}\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare DeepImmuno input from netMHCpan output.")
    p.add_argument("--input-csv", required=True, help="Input netMHCpan filtered CSV")
    p.add_argument("--output-input", required=True, help="Output DeepImmuno input CSV")
    p.add_argument("--output-map", required=True, help="Output metadata map CSV")
    p.add_argument("--pair-id", required=True, help="Pair ID")
    p.add_argument("--lengths", default="9,10", help="Comma-separated supported peptide lengths")
    p.add_argument("--log", required=True, help="Log file path")
    return p.parse_args()


def run_from_snakemake() -> None:
    pair_id = str(getattr(snakemake.wildcards, "pair_id", "unknown"))
    prepare_inputs(
        in_csv=Path(str(snakemake.input.netmhcpan_filtered)),
        out_input_csv=Path(str(snakemake.output.input_csv)),
        out_map_csv=Path(str(snakemake.output.map_csv)),
        pair_id=pair_id,
        lengths=parse_lengths(snakemake.params.lengths),
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.input_csv)).with_suffix(".log"),
    )


def main() -> None:
    args = parse_args()
    prepare_inputs(
        in_csv=Path(args.input_csv),
        out_input_csv=Path(args.output_input),
        out_map_csv=Path(args.output_map),
        pair_id=args.pair_id,
        lengths=parse_lengths(args.lengths),
        log_path=Path(args.log),
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
