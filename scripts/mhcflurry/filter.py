#!/usr/bin/env python3
"""Filter MHCflurry predictions by percentile-like rank threshold."""

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


def pick_rank_column(fieldnames: list[str]) -> str:
    lowered = {name.lower(): name for name in fieldnames}
    preferred = [
        "presentation_percentile",
        "presentation_percentile_rank",
        "percentile_rank",
        "affinity_percentile",
    ]
    for col in preferred:
        if col in lowered:
            return lowered[col]
    for field in fieldnames:
        if "percentile" in field.lower():
            return field
    raise ValueError(f"No percentile/rank-like column found in fields: {fieldnames}")


def filter_rows(rows: list[dict[str, str]], rank_column: str, rank_threshold: float) -> tuple[list[dict[str, str]], int]:
    kept: list[dict[str, str]] = []
    invalid = 0
    for row in rows:
        rank_value = to_float(row.get(rank_column))
        if rank_value is None:
            invalid += 1
            continue
        if rank_value <= rank_threshold:
            kept.append(row)
    return kept, invalid


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_log(
    path: Path,
    rank_column: str,
    rank_threshold: float,
    total: int,
    kept: int,
    invalid: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as log:
        log.write("MHCflurry filter summary\n")
        log.write(f"rank_column: {rank_column}\n")
        log.write(f"rank_threshold: {rank_threshold}\n")
        log.write(f"total_rows: {total}\n")
        log.write(f"kept_rows: {kept}\n")
        log.write(f"filtered_out_rows: {total - kept}\n")
        log.write(f"invalid_rank_rows: {invalid}\n")


def run_pipeline(
    in_csv: Path,
    out_csv: Path,
    log_path: Path,
    rank_threshold: float,
) -> None:
    fieldnames, rows = read_rows(in_csv)
    rank_column = pick_rank_column(fieldnames)
    kept, invalid = filter_rows(rows, rank_column, rank_threshold)
    write_rows(out_csv, fieldnames, kept)
    write_log(log_path, rank_column, rank_threshold, len(rows), len(kept), invalid)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Filter MHCflurry predictions by rank threshold.")
    p.add_argument("--input", required=True, help="Input all.csv")
    p.add_argument("--output", required=True, help="Output filtered.csv")
    p.add_argument("--log", required=True, help="Log file path")
    p.add_argument("--rank-threshold", type=float, required=True, help="Maximum rank/percentile")
    return p.parse_args()


def run_from_snakemake() -> None:
    run_pipeline(
        in_csv=Path(str(snakemake.input.all_csv)),
        out_csv=Path(str(snakemake.output.filtered)),
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.filtered)).with_suffix(".log"),
        rank_threshold=float(snakemake.params.rank_threshold),
    )


def main() -> None:
    args = parse_args()
    run_pipeline(
        in_csv=Path(args.input),
        out_csv=Path(args.output),
        log_path=Path(args.log),
        rank_threshold=args.rank_threshold,
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
