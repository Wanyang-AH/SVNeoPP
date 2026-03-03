#!/usr/bin/env python3
"""Generate FragPipe manifest from proteomics rows in samples.tsv."""

import csv
from pathlib import Path


def _int_or_default(value, default):
    try:
        return int(str(value).strip())
    except Exception:
        return default


def _status_rank(status):
    text = (status or "").strip().lower()
    if text == "tumor":
        return 0
    if text == "normal":
        return 1
    return 9


def _container_spec_path(raw_path, mount_prefix):
    raw = str(raw_path).strip()
    if raw.startswith("/"):
        return raw
    return f"{mount_prefix.rstrip('/')}/{raw.lstrip('./')}"


def build_manifest(samples_tsv, pair_id, out_manifest, log_path, mount_prefix):
    rows = []
    with samples_tsv.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if (row.get("datatype") or "").strip().lower() != "prot":
                continue
            if (row.get("pair_id") or "").strip() != pair_id:
                continue
            rows.append(row)

    rows.sort(
        key=lambda r: (
            _status_rank(r.get("status")),
            _int_or_default(r.get("replicate"), 99999),
            (r.get("sample_id") or "").strip(),
        )
    )

    if not rows:
        raise ValueError(f"No proteomics rows found for pair_id={pair_id} in {samples_tsv}")

    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    missing_files = 0
    with out_manifest.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["spec_file", "experiment", "bioreplicate", "condition"])

        for idx, row in enumerate(rows, start=1):
            raw_file = (row.get("r1") or "").strip()
            if not raw_file:
                raise ValueError(f"Empty r1 for pair_id={pair_id}, sample={row.get('sample_id', 'NA')}")
            if not Path(raw_file).exists():
                missing_files += 1

            sample_id = (row.get("sample_id") or f"{pair_id}_prot_{idx}").strip()
            replicate = str(_int_or_default(row.get("replicate"), idx))
            condition = (row.get("status") or "unknown").strip().lower()
            spec_file = _container_spec_path(raw_file, mount_prefix)

            writer.writerow([spec_file, sample_id, replicate, condition])

    with log_path.open("w") as log_handle:
        log_handle.write("FragPipe manifest preparation summary\n")
        log_handle.write(f"samples_tsv: {samples_tsv}\n")
        log_handle.write(f"pair_id: {pair_id}\n")
        log_handle.write(f"output_manifest: {out_manifest}\n")
        log_handle.write(f"rows_written: {len(rows)}\n")
        log_handle.write(f"missing_raw_files: {missing_files}\n")
        log_handle.write(f"mount_prefix: {mount_prefix}\n")

    if missing_files:
        raise ValueError(f"{missing_files} RAW file(s) missing for pair_id={pair_id}. See {log_path}.")


def run_from_snakemake():
    samples_tsv = Path(str(snakemake.input.samples_tsv))
    out_manifest = Path(str(snakemake.output.manifest))
    log_path = Path(str(snakemake.log[0]))
    pair_id = str(snakemake.wildcards.pair_id)
    mount_prefix = str(snakemake.params.mount_prefix)
    build_manifest(samples_tsv, pair_id, out_manifest, log_path, mount_prefix)


if __name__ == "__main__":
    run_from_snakemake()
