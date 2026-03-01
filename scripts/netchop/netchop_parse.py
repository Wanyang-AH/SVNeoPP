#!/usr/bin/env python3
"""Parse NetChop raw output and compute peptide-level cleavage metrics."""

import argparse
import csv
import re
from pathlib import Path

ROW_RE = re.compile(
    r"^\s*(\d+)\s+([A-Za-z])\s+([A-Za-z\.])\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
)

METRIC_COLUMNS = [
    "pair_id",
    "ident",
    "sv_id",
    "gene_name",
    "peptide_source",
    "sv_type",
    "sample",
    "peptide_sequence",
    "peptide_length",
    "N_score",
    "C_score",
    "max_internal_score",
    "internal_cleavage_count",
]


def parse_fasta_entries(fasta_path: Path) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    header = ""
    with fasta_path.open("r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]
                continue

            parts = header.split("|")
            entry = {
                "sv_id": parts[0] if len(parts) > 0 else "NA",
                "gene_name": parts[1] if len(parts) > 1 else "NA",
                "peptide_source": parts[2] if len(parts) > 2 else "NA",
                "sv_type": parts[3] if len(parts) > 3 else "NA",
                "sample": parts[4] if len(parts) > 4 else "NA",
                "ident": parts[5] if len(parts) > 5 else f"idx{len(entries) + 1}",
                "peptide_sequence": line,
                "peptide_length": len(line),
            }
            entries.append(entry)
            header = ""
    return entries


def compute_metrics(block: list[tuple[int, float, bool]]) -> dict[str, float | int]:
    if not block:
        return {
            "N_score": 0.0,
            "C_score": 0.0,
            "max_internal_score": 0.0,
            "internal_cleavage_count": 0,
        }

    n_score = block[0][1]
    c_score = block[-1][1]
    internal = block[1:-1]

    if internal:
        max_internal = max(x[1] for x in internal)
        internal_count = sum(1 for x in internal if x[2])
    else:
        max_internal = 0.0
        internal_count = 0

    return {
        "N_score": n_score,
        "C_score": c_score,
        "max_internal_score": max_internal,
        "internal_cleavage_count": internal_count,
    }


def parse_netchop_blocks(raw_path: Path) -> list[dict[str, float | int]]:
    blocks: list[dict[str, float | int]] = []
    current: list[tuple[int, float, bool]] = []
    last_pos = 0

    def flush() -> None:
        nonlocal current, last_pos
        if current:
            blocks.append(compute_metrics(current))
            current = []
            last_pos = 0

    with raw_path.open("r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith("NetChop"):
                continue
            if set(line) == {"-"} or line.startswith("----"):
                flush()
                continue

            m = ROW_RE.match(line)
            if not m:
                continue

            pos = int(m.group(1))
            marker = m.group(3)
            score = float(m.group(4))
            is_cleavage = marker.upper() == "S"

            if pos == 1 and current and last_pos >= 1:
                flush()
            current.append((pos, score, is_cleavage))
            last_pos = pos

    flush()
    return blocks


def write_metrics(
    fasta_entries: list[dict[str, str]],
    netchop_metrics: list[dict[str, float | int]],
    out_csv: Path,
    pair_id: str,
    log_path: Path,
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    matched = min(len(fasta_entries), len(netchop_metrics))
    rows: list[dict[str, str | float | int]] = []
    for idx in range(matched):
        merged = {"pair_id": pair_id, **fasta_entries[idx], **netchop_metrics[idx]}
        rows.append(merged)

    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METRIC_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    with log_path.open("w") as log:
        log.write("NetChop parse summary\n")
        log.write(f"pair_id: {pair_id}\n")
        log.write(f"fasta_entries: {len(fasta_entries)}\n")
        log.write(f"netchop_blocks: {len(netchop_metrics)}\n")
        log.write(f"matched_rows: {matched}\n")
        if len(fasta_entries) != len(netchop_metrics):
            log.write("warning: FASTA and NetChop counts differ, truncated to min length.\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Parse NetChop raw output to metrics CSV.")
    p.add_argument("--fasta", required=True, help="Input FASTA used for netChop")
    p.add_argument("--raw", required=True, help="Raw netChop output")
    p.add_argument("--output", required=True, help="Output metrics CSV")
    p.add_argument("--pair-id", required=True, help="Pair/sample ID")
    p.add_argument("--log", required=True, help="Log file")
    return p.parse_args()


def run_pipeline(fasta: Path, raw: Path, out_csv: Path, pair_id: str, log: Path) -> None:
    fasta_entries = parse_fasta_entries(fasta)
    metrics = parse_netchop_blocks(raw)
    write_metrics(fasta_entries, metrics, out_csv, pair_id, log)


def run_from_snakemake() -> None:
    fasta = Path(str(snakemake.input.fasta))
    raw = Path(str(snakemake.input.raw))
    out_csv = Path(str(snakemake.output.metrics))
    pair_id = str(getattr(snakemake.wildcards, "pair_id", "unknown"))
    log = Path(str(snakemake.log[0])) if snakemake.log else out_csv.with_suffix(".log")
    run_pipeline(fasta, raw, out_csv, pair_id, log)


def main() -> None:
    args = parse_args()
    run_pipeline(
        fasta=Path(args.fasta),
        raw=Path(args.raw),
        out_csv=Path(args.output),
        pair_id=args.pair_id,
        log=Path(args.log),
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
