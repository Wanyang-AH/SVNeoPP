#!/usr/bin/env python3
"""Run DeepImmuno-CNN in multiple mode for one pair."""

import argparse
import shlex
import subprocess
from pathlib import Path


def has_content(path: Path) -> bool:
    with path.open("r") as handle:
        for line in handle:
            if line.strip():
                return True
    return False


def run_deepimmuno(
    input_csv: Path,
    output_tsv: Path,
    deepimmuno_home: Path,
    python_bin: str,
    extra: str,
    log_path: Path,
) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    deep_script = deepimmuno_home / "deepimmuno-cnn.py"
    if not deep_script.exists():
        raise FileNotFoundError(f"DeepImmuno script not found: {deep_script}")

    if not has_content(input_csv):
        with output_tsv.open("w") as out:
            out.write("peptide\tHLA\timmunogenicity\n")
        with log_path.open("w") as log:
            log.write("DeepImmuno run skipped: no eligible 9/10-mer entries.\n")
            log.write(f"input_csv: {input_csv}\n")
            log.write(f"output_tsv: {output_tsv}\n")
        return

    out_dir = output_tsv.parent.resolve()
    raw_result = out_dir / "deepimmuno-cnn-result.txt"
    if raw_result.exists():
        raw_result.unlink()

    cmd = [
        python_bin,
        str(deep_script.resolve()),
        "--mode",
        "multiple",
        "--intdir",
        str(input_csv.resolve()),
        "--outdir",
        str(out_dir),
    ]
    if extra.strip():
        cmd.extend(shlex.split(extra))

    with log_path.open("w") as log:
        log.write(f"DeepImmuno command: {' '.join(cmd)}\n")
        log.write(f"DeepImmuno home: {deepimmuno_home}\n")
        subprocess.run(
            cmd,
            cwd=str(deepimmuno_home),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )

    if not raw_result.exists():
        raise FileNotFoundError(f"Expected DeepImmuno output not found: {raw_result}")

    if output_tsv.exists():
        output_tsv.unlink()
    raw_result.replace(output_tsv)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run DeepImmuno-CNN on prepared peptide/HLA input.")
    p.add_argument("--input", required=True, help="Input CSV with peptide,HLA rows (no header)")
    p.add_argument("--output", required=True, help="Output TSV path")
    p.add_argument("--deepimmuno-home", required=True, help="Path to DeepImmuno directory")
    p.add_argument("--python-bin", default="python", help="Python executable used to run DeepImmuno")
    p.add_argument("--extra", default="", help="Extra CLI args passed to deepimmuno-cnn.py")
    p.add_argument("--log", required=True, help="Log file path")
    return p.parse_args()


def run_from_snakemake() -> None:
    run_deepimmuno(
        input_csv=Path(str(snakemake.input.input_csv)),
        output_tsv=Path(str(snakemake.output.raw_tsv)),
        deepimmuno_home=Path(str(snakemake.params.deepimmuno_home)),
        python_bin=str(snakemake.params.python_bin),
        extra=str(snakemake.params.extra),
        log_path=Path(str(snakemake.log[0])) if snakemake.log else Path(str(snakemake.output.raw_tsv)).with_suffix(".log"),
    )


def main() -> None:
    args = parse_args()
    run_deepimmuno(
        input_csv=Path(args.input),
        output_tsv=Path(args.output),
        deepimmuno_home=Path(args.deepimmuno_home),
        python_bin=args.python_bin,
        extra=args.extra,
        log_path=Path(args.log),
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
