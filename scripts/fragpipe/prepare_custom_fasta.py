#!/usr/bin/env python3
"""Combine neoantigen FASTA and reference FASTA for FragPipe Philosopher database step."""

from pathlib import Path


def count_fasta_records(path):
    count = 0
    with path.open("r") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def copy_fasta_content(source, target_handle):
    with source.open("r") as in_handle:
        for line in in_handle:
            target_handle.write(line)
    if source.stat().st_size > 0:
        with source.open("rb") as in_handle:
            in_handle.seek(-1, 2)
            if in_handle.read(1) != b"\n":
                target_handle.write("\n")


def build_custom_fasta(neo_fasta, reference_fasta, output_fasta, log_path):
    if not neo_fasta.exists():
        raise FileNotFoundError(f"Neo FASTA not found: {neo_fasta}")
    if not reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    neo_count = count_fasta_records(neo_fasta)
    reference_count = count_fasta_records(reference_fasta)

    with output_fasta.open("w") as out_handle:
        copy_fasta_content(neo_fasta, out_handle)
        copy_fasta_content(reference_fasta, out_handle)

    out_count = count_fasta_records(output_fasta)

    with log_path.open("w") as log_handle:
        log_handle.write("FragPipe custom FASTA preparation summary\n")
        log_handle.write(f"neo_fasta: {neo_fasta}\n")
        log_handle.write(f"reference_fasta: {reference_fasta}\n")
        log_handle.write(f"output_fasta: {output_fasta}\n")
        log_handle.write(f"neo_records: {neo_count}\n")
        log_handle.write(f"reference_records: {reference_count}\n")
        log_handle.write(f"output_records: {out_count}\n")

    if out_count < neo_count + reference_count:
        raise ValueError(
            "Output FASTA record count is lower than input sum. "
            f"neo={neo_count}, reference={reference_count}, output={out_count}"
        )


def run_from_snakemake():
    neo_fasta = Path(str(snakemake.input.neo_fasta))
    reference_fasta = Path(str(snakemake.input.reference_fasta))
    output_fasta = Path(str(snakemake.output.custom_fasta))
    log_path = Path(str(snakemake.log[0]))
    build_custom_fasta(neo_fasta, reference_fasta, output_fasta, log_path)


if __name__ == "__main__":
    run_from_snakemake()
