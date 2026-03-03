#!/usr/bin/env python3
"""Patch FragPipe workflow database path for per-pair generated decoy FASTA."""

from pathlib import Path


def patch_workflow(template_path, output_path, db_path_in_container, log_path):
    if not template_path.exists():
        raise FileNotFoundError(f"Workflow template not found: {template_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    replaced = False
    out_lines = []
    with template_path.open("r") as handle:
        for line in handle:
            if line.startswith("database.db-path="):
                out_lines.append(f"database.db-path={db_path_in_container}\n")
                replaced = True
            else:
                out_lines.append(line)

    if not replaced:
        out_lines.insert(0, f"database.db-path={db_path_in_container}\n")

    with output_path.open("w") as out_handle:
        out_handle.writelines(out_lines)

    with log_path.open("w") as log_handle:
        log_handle.write("FragPipe workflow patch summary\n")
        log_handle.write(f"template: {template_path}\n")
        log_handle.write(f"output: {output_path}\n")
        log_handle.write(f"database_db_path: {db_path_in_container}\n")
        log_handle.write(f"database_line_replaced: {replaced}\n")


def run_from_snakemake():
    template_path = Path(str(snakemake.input.template))
    output_path = Path(str(snakemake.output.workflow))
    log_path = Path(str(snakemake.log[0]))
    db_path_in_container = str(snakemake.params.db_path_in_container)
    patch_workflow(template_path, output_path, db_path_in_container, log_path)


if __name__ == "__main__":
    run_from_snakemake()
