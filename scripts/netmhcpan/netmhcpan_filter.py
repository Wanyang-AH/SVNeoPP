#!/usr/bin/env python3
"""Filter netMHCpan predictions by EL and BA rank thresholds."""

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


def filter_rows(
    rows: list[dict[str, str]],
    rank_el_threshold: float,
    rank_ba_threshold: float,
) -> tuple[list[dict[str, str]], int]:
    kept: list[dict[str, str]] = []
    invalid = 0
    for row in rows:
        el_rank = to_float(row.get("el_rank"))
        ba_rank = to_float(row.get("ba_rank"))
        if el_rank is None or ba_rank is None:
            invalid += 1
            continue
        if el_rank <= rank_el_threshold and ba_rank <= rank_ba_threshold:
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
    total: int,
    kept: int,
    invalid: int,
    rank_el_threshold: float,
    rank_ba_threshold: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as log:
        log.write("netMHCpan filter summary\n")
        log.write(f"total_rows: {total}\n")
        log.write(f"kept_rows: {kept}\n")
        log.write(f"filtered_out_rows: {total - kept}\n")
        log.write(f"invalid_rank_rows: {invalid}\n")
        log.write(f"rank_el_threshold: {rank_el_threshold}\n")
        log.write(f"rank_ba_threshold: {rank_ba_threshold}\n")


def run_pipeline(
    in_csv: Path,
    out_csv: Path,
    log_path: Path,
    rank_el_threshold: float,
    rank_ba_threshold: float,
) -> None:
    fieldnames, rows = read_rows(in_csv)
    kept, invalid = filter_rows(rows, rank_el_threshold, rank_ba_threshold)
    write_rows(out_csv, fieldnames, kept)
    write_log(log_path, len(rows), len(kept), invalid, rank_el_threshold, rank_ba_threshold)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Filter netMHCpan predictions using EL/BA rank thresholds.")
    p.add_argument("--input", required=True, help="Input all.csv")
    p.add_argument("--output", required=True, help="Output filtered.csv")
    p.add_argument("--log", required=True, help="Log file path")
    p.add_argument("--rank-el-threshold", type=float, required=True, help="Maximum %Rank_EL")
    p.add_argument("--rank-ba-threshold", type=float, required=True, help="Maximum %Rank_BA")
    return p.parse_args()


def run_from_snakemake() -> None:
    run_pipeline(
        in_csv=Path(str(snakemake.input.all_csv)),
        out_csv=Path(str(snakemake.output.filtered)),
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.filtered)).with_suffix(".log"),
        rank_el_threshold=float(snakemake.params.rank_el_threshold),
        rank_ba_threshold=float(snakemake.params.rank_ba_threshold),
    )


def main() -> None:
    args = parse_args()
    run_pipeline(
        in_csv=Path(args.input),
        out_csv=Path(args.output),
        log_path=Path(args.log),
        rank_el_threshold=args.rank_el_threshold,
        rank_ba_threshold=args.rank_ba_threshold,
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
