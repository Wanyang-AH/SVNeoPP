from pathlib import Path
import sys

import pandas as pd


def read_csv(csv_path: Path) -> pd.DataFrame:
    """Load netchop CSV and validate required columns."""
    if not csv_path.exists():
        sys.exit(f"Missing input file: {csv_path}")

    df = pd.read_csv(csv_path)
    required = ["Ident", "peptide_sequence", "peptide_length", "sample", "Gene"]
    missing_cols = [col for col in required if col not in df.columns]
    if missing_cols:
        sys.exit(f"CSV is missing required columns: {', '.join(missing_cols)}")

    return df


def write_fastas(df: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    n_missing_seq = df["peptide_sequence"].isna().sum()
    if n_missing_seq:
        print(f"Skipping {n_missing_seq} rows with missing peptide_sequence", file=sys.stderr)
        df = df[df["peptide_sequence"].notna()].copy()

    df["Gene"] = df["Gene"].fillna("NA")

    # 仅保留长度 8-11 的肽段，按样本+长度拆分
    valid_lengths = {8, 9, 10, 11}
    df = df[df["peptide_length"].isin(valid_lengths)].copy()

    for (sample, length), group in df.groupby(["sample", "peptide_length"]):
        fasta_lines = []
        for idx, row in group.iterrows():
            header = (
                f">{idx + 1}|{row['sample']}|{row['Gene']}|"
                f"{row['peptide_sequence']}|{row['peptide_length']}"
            )
            fasta_lines.append(header)
            fasta_lines.append(str(row["peptide_sequence"]))

        output_path = output_dir / f"{sample}_len{length}.fasta"
        output_path.write_text("\n".join(fasta_lines) + "\n")
        print(f"Wrote {len(group)} entries to {output_path}")


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    csv_path = repo_root / "1120_netchop" / "netchop.csv"
    output_dir = repo_root / "1121mhc"

    df = read_csv(csv_path)
    write_fastas(df, output_dir)


if __name__ == "__main__":
    main()
