#!/usr/bin/env python3
"""Filter DeepImmuno predictions by immunogenicity threshold."""

import argparse
import csv
from pathlib import Path


def to_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return fieldnames, rows


def filter_rows(rows: list[dict[str, str]], threshold: float) -> tuple[list[dict[str, str]], int]:
    kept: list[dict[str, str]] = []
    invalid = 0
    for row in rows:
        score = to_float(row.get("deepimmuno_score"))
        if score is None:
            invalid += 1
            continue
        if score >= threshold:
            kept.append(row)
    return kept, invalid


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_log(path: Path, threshold: float, total: int, kept: int, invalid: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as log:
        log.write("DeepImmuno filter summary\n")
        log.write(f"threshold: {threshold}\n")
        log.write(f"total_rows: {total}\n")
        log.write(f"kept_rows: {kept}\n")
        log.write(f"filtered_out_rows: {total - kept}\n")
        log.write(f"invalid_score_rows: {invalid}\n")


def run_pipeline(in_csv: Path, out_csv: Path, log_path: Path, threshold: float) -> None:
    fieldnames, rows = read_rows(in_csv)
    if "deepimmuno_score" not in fieldnames:
        raise ValueError(f"Missing required column deepimmuno_score in {in_csv}")
    kept, invalid = filter_rows(rows, threshold)
    write_rows(out_csv, fieldnames, kept)
    write_log(log_path, threshold, len(rows), len(kept), invalid)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Filter DeepImmuno predictions using threshold.")
    p.add_argument("--input", required=True, help="Input all.csv")
    p.add_argument("--output", required=True, help="Output filtered.csv")
    p.add_argument("--log", required=True, help="Log file path")
    p.add_argument("--threshold", required=True, type=float, help="Minimum deepimmuno_score to keep")
    return p.parse_args()


def run_from_snakemake() -> None:
    run_pipeline(
        in_csv=Path(str(snakemake.input.all_csv)),
        out_csv=Path(str(snakemake.output.filtered)),
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.filtered)).with_suffix(".log"),
        threshold=float(snakemake.params.threshold),
    )


def main() -> None:
    args = parse_args()
    run_pipeline(
        in_csv=Path(args.input),
        out_csv=Path(args.output),
        log_path=Path(args.log),
        threshold=args.threshold,
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
