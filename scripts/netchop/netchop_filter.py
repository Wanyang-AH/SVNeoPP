#!/usr/bin/env python3
"""Filter NetChop metrics to retain peptide candidates."""

import argparse
import csv
from pathlib import Path


def to_float(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def to_int(value: str | None, default: int = 0) -> int:
    try:
        return int(float(value)) if value not in (None, "") else default
    except ValueError:
        return default


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return fieldnames, rows


def filter_rows(
    rows: list[dict[str, str]],
    n_term_min: float,
    c_term_min: float,
    internal_max: float,
    internal_count_max: int,
) -> list[dict[str, str]]:
    kept: list[dict[str, str]] = []
    for row in rows:
        n_score = to_float(row.get("N_score"))
        c_score = to_float(row.get("C_score"))
        max_internal = to_float(row.get("max_internal_score"))
        internal_count = to_int(row.get("internal_cleavage_count"))

        if n_score < n_term_min:
            continue
        if c_score < c_term_min:
            continue
        if max_internal > internal_max:
            continue
        if internal_count > internal_count_max:
            continue
        kept.append(row)
    return kept


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_log(
    path: Path,
    total: int,
    kept: int,
    n_term_min: float,
    c_term_min: float,
    internal_max: float,
    internal_count_max: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as log:
        log.write("NetChop filter summary\n")
        log.write(f"total_rows: {total}\n")
        log.write(f"kept_rows: {kept}\n")
        log.write(f"dropped_rows: {total - kept}\n")
        log.write(f"n_term_min: {n_term_min}\n")
        log.write(f"c_term_min: {c_term_min}\n")
        log.write(f"internal_max: {internal_max}\n")
        log.write(f"internal_count_max: {internal_count_max}\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Filter NetChop metrics table.")
    p.add_argument("--input", required=True, help="Input metrics CSV")
    p.add_argument("--output", required=True, help="Output filtered CSV")
    p.add_argument("--log", required=True, help="Output log path")
    p.add_argument("--n-term-min", type=float, required=True, help="Minimum N_score")
    p.add_argument("--c-term-min", type=float, required=True, help="Minimum C_score")
    p.add_argument("--internal-max", type=float, required=True, help="Maximum internal score")
    p.add_argument("--internal-count-max", type=int, required=True, help="Maximum internal cleavage count")
    return p.parse_args()


def run_pipeline(
    in_csv: Path,
    out_csv: Path,
    log_path: Path,
    n_term_min: float,
    c_term_min: float,
    internal_max: float,
    internal_count_max: int,
) -> None:
    fieldnames, rows = read_rows(in_csv)
    kept = filter_rows(rows, n_term_min, c_term_min, internal_max, internal_count_max)
    write_rows(out_csv, fieldnames, kept)
    write_log(log_path, len(rows), len(kept), n_term_min, c_term_min, internal_max, internal_count_max)


def run_from_snakemake() -> None:
    run_pipeline(
        in_csv=Path(str(snakemake.input.metrics)),
        out_csv=Path(str(snakemake.output.filtered)),
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.filtered)).with_suffix(".log"),
        n_term_min=float(snakemake.params.n_term_min),
        c_term_min=float(snakemake.params.c_term_min),
        internal_max=float(snakemake.params.internal_max),
        internal_count_max=int(snakemake.params.internal_count_max),
    )


def main() -> None:
    args = parse_args()
    run_pipeline(
        in_csv=Path(args.input),
        out_csv=Path(args.output),
        log_path=Path(args.log),
        n_term_min=args.n_term_min,
        c_term_min=args.c_term_min,
        internal_max=args.internal_max,
        internal_count_max=args.internal_count_max,
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
