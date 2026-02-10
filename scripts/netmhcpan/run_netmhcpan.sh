#!/bin/bash

# Run netMHCpan on per-sample, per-length FASTA files in parallel (up to 8 concurrent jobs).
# Expects input FASTA files named <sample>_len<8|9|10|11>.fasta under INPUT_DIR.
# Output files are written to OUTPUT_DIR as <sample>_len<len>_netmhcpan.txt.

set -euo pipefail

# ================= 配置区域 =================
INPUT_DIR="/sibpt/anwy/lungfile/1121mhc"
OUTPUT_DIR="${INPUT_DIR}/netMHCpan"
MAX_JOBS=8

mkdir -p "$OUTPUT_DIR"

# 样本对应 HLA
declare -A SAMPLES
# netMHCpan 4.x 接受不带 "*" 的等位基因格式，如 HLA-A02:01
SAMPLES["l041"]="HLA-A02:01,HLA-A24:02,HLA-B46:01,HLA-B15:01,HLA-C01:02,HLA-C08:01"
SAMPLES["l048"]="HLA-A02:07,HLA-A23:01,HLA-B45:01,HLA-B46:01,HLA-C01:02,HLA-C06:02"
SAMPLES["l052"]="HLA-A33:03,HLA-A24:02,HLA-B58:01,HLA-B40:01,HLA-C03:02,HLA-C03:04"
SAMPLES["l056"]="HLA-A02:07,HLA-B46:01,HLA-C01:02"

lengths=(8 9 10 11)

run_cmd() {
  local sample="$1"
  local len="$2"
  local alleles="$3"
  local infile="${INPUT_DIR}/${sample}_len${len}.fasta"
  local outfile="${OUTPUT_DIR}/${sample}_len${len}_netmhcpan.txt"

  if [[ ! -f "$infile" ]]; then
    echo "[WARN] Missing input: $infile" >&2
    return
  fi

  echo "[RUN] sample=${sample} len=${len} -> $outfile"
  netMHCpan -f "$infile" -l "$len" -a "$alleles" -s -BA > "$outfile"
}

echo "Starting netMHCpan runs (max ${MAX_JOBS} concurrent)..."

for sample in "${!SAMPLES[@]}"; do
  alleles="${SAMPLES[$sample]}"
  for len in "${lengths[@]}"; do
    # 控制并发数
    while [[ $(jobs -r | wc -l) -ge $MAX_JOBS ]]; do
      sleep 0.5
    done
    run_cmd "$sample" "$len" "$alleles" &
  done
done

wait
echo "All netMHCpan jobs completed. Outputs in $OUTPUT_DIR"
