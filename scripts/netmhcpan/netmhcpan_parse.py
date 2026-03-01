#!/usr/bin/env python3
"""Parse netMHCpan xls-style output into a normalized long-format CSV."""

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


def read_map(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        out: dict[str, dict[str, str]] = {}
        for row in reader:
            peptide = (row.get("peptide_sequence") or "").strip()
            if peptide and peptide not in out:
                out[peptide] = row
    return fieldnames, out


def parse_xls_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]

    if len(lines) < 3:
        raise ValueError(f"netMHCpan xls file is too short: {path}")

    header_alleles = [x.strip() for x in lines[0].split("\t")]
    allele_blocks: list[tuple[int, str]] = []
    for idx, cell in enumerate(header_alleles):
        if cell.startswith("HLA-"):
            allele_blocks.append((idx, cell))
    if not allele_blocks:
        raise ValueError(f"No HLA allele blocks found in xls header: {path}")

    header_names = [x.strip() for x in lines[1].split("\t")]
    data_rows = lines[2:]

    records: list[dict[str, str]] = []
    for raw in data_rows:
        cols = [x.strip() for x in raw.split("\t")]
        if len(cols) < 2:
            continue
        peptide = cols[1] if len(cols) > 1 else ""
        if not peptide:
            continue

        ave = cols[-2] if len(cols) >= 2 else ""
        nb = cols[-1] if len(cols) >= 1 else ""
        row_id = cols[2] if len(cols) > 2 else ""
        meas = cols[3] if len(cols) > 3 else ""
        pos = cols[0] if len(cols) > 0 else ""

        for start_idx, allele in allele_blocks:
            end_idx = start_idx + 5
            if end_idx >= len(cols):
                continue
            rec = {
                "pos": pos,
                "peptide": peptide,
                "id": row_id,
                "meas": meas,
                "allele": allele,
                "core": cols[start_idx],
                "icore": cols[start_idx + 1],
                "el_score": cols[start_idx + 2],
                "el_rank": cols[start_idx + 3],
                "ba_score": cols[start_idx + 4],
                "ba_rank": cols[start_idx + 5],
                "ave": ave,
                "nb": nb,
            }
            records.append(rec)

    return header_names, records


def write_all_csv(
    out_csv: Path,
    pair_id: str,
    parsed_rows: list[dict[str, str]],
    map_fieldnames: list[str],
    peptide_map: dict[str, dict[str, str]],
) -> int:
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    map_fields = [x for x in map_fieldnames if x != "peptide_sequence"]
    base_fields = [
        "pair_id",
        "peptide_sequence",
        "peptide_length",
        "allele",
        "pos",
        "id",
        "meas",
        "core",
        "icore",
        "el_score",
        "el_rank",
        "ba_score",
        "ba_rank",
        "ave",
        "nb",
    ]
    out_fields = base_fields + [x for x in map_fields if x not in base_fields]

    n_with_map = 0
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields)
        writer.writeheader()
        for row in parsed_rows:
            peptide = row["peptide"]
            map_row = peptide_map.get(peptide, {})
            if map_row:
                n_with_map += 1

            out_row = {
                "pair_id": pair_id,
                "peptide_sequence": peptide,
                "peptide_length": len(peptide),
                "allele": row["allele"],
                "pos": row["pos"],
                "id": row["id"],
                "meas": row["meas"],
                "core": row["core"],
                "icore": row["icore"],
                "el_score": row["el_score"],
                "el_rank": row["el_rank"],
                "ba_score": row["ba_score"],
                "ba_rank": row["ba_rank"],
                "ave": row["ave"],
                "nb": row["nb"],
            }
            for key in map_fields:
                out_row[key] = map_row.get(key, "")
            writer.writerow(out_row)

    return n_with_map


def write_log(
    path: Path,
    pair_id: str,
    xls_path: Path,
    map_path: Path,
    n_rows: int,
    n_with_map: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as log:
        log.write("netMHCpan parse summary\n")
        log.write(f"pair_id: {pair_id}\n")
        log.write(f"xls_input: {xls_path}\n")
        log.write(f"map_input: {map_path}\n")
        log.write(f"parsed_rows: {n_rows}\n")
        log.write(f"rows_with_metadata: {n_with_map}\n")


def run_pipeline(
    xls_path: Path,
    map_path: Path,
    out_csv: Path,
    pair_id: str,
    log_path: Path,
) -> None:
    _, parsed_rows = parse_xls_table(xls_path)
    map_fieldnames, peptide_map = read_map(map_path)
    n_with_map = write_all_csv(out_csv, pair_id, parsed_rows, map_fieldnames, peptide_map)
    write_log(log_path, pair_id, xls_path, map_path, len(parsed_rows), n_with_map)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Parse netMHCpan xls output to long-format CSV.")
    p.add_argument("--input-xls", required=True, help="netMHCpan xls-style output file")
    p.add_argument("--input-map", required=True, help="Peptide metadata map CSV")
    p.add_argument("--output", required=True, help="Output CSV")
    p.add_argument("--pair-id", required=True, help="Pair ID")
    p.add_argument("--log", required=True, help="Log file")
    return p.parse_args()


def run_from_snakemake() -> None:
    pair_id = str(getattr(snakemake.wildcards, "pair_id", "unknown"))
    run_pipeline(
        xls_path=Path(str(snakemake.input.xls)),
        map_path=Path(str(snakemake.input.map_csv)),
        out_csv=Path(str(snakemake.output.all_csv)),
        pair_id=pair_id,
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.all_csv)).with_suffix(".log"),
    )


def main() -> None:
    args = parse_args()
    run_pipeline(
        xls_path=Path(args.input_xls),
        map_path=Path(args.input_map),
        out_csv=Path(args.output),
        pair_id=args.pair_id,
        log_path=Path(args.log),
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
