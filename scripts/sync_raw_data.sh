#!/usr/bin/env bash
set -euo pipefail

REMOTE_HOST="anwy@172.16.50.2"
REMOTE_BASE="/mnt/data/lumanman/result_20250603_all"
PORT="22"
PARALLEL=2

WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v sshpass >/dev/null 2>&1; then
  echo "ERROR: sshpass is not installed. Install it or switch to SSH key auth."
  exit 1
fi

read -s -p "SSH password for ${REMOTE_HOST}: " SSHPASS
echo

sync_one() {
  local remote_path="$1"
  local local_path="$2"
  mkdir -p "$(dirname "${local_path}")"
  sshpass -p "${SSHPASS}" rsync -av --progress --partial --append-verify \
    -e "ssh -p ${PORT}" \
    "${REMOTE_HOST}:${remote_path}" "${local_path}"
}

# Export for xargs workers
export SSHPASS REMOTE_HOST PORT
export -f sync_one

TASKS=()

# WGS (DNA)
for id in 041 048 052 056; do
  for status in normal tumor; do
    for read in R1 R2; do
      remote="${REMOTE_BASE}/L${id}/dna/00_merged/${status}/L${id}_dna_${status}_${read}.fastq.gz"
      local="${WORKDIR}/data/wgs/l${id}/l${id}_dna_${status}_${read}.fastq.gz"
      TASKS+=("${remote}"$'\t'"${local}")
    done
  done
done

# RNA
for id in 041 048 052 056; do
  for status in normal tumor; do
    for read in R1 R2; do
      remote="${REMOTE_BASE}/L${id}/rna/00_merged/${status}/L${id}_rna_${status}_${read}.fastq.gz"
      local="${WORKDIR}/data/rna/l${id}/l${id}_rna_${status}_${read}.fastq.gz"
      TASKS+=("${remote}"$'\t'"${local}")
    done
  done
done

# Proteomics (RAW) — rename L041 normal P -> N to match samples.tsv
while IFS= read -r remote; do
  [ -z "${remote}" ] && continue
  if [[ "${remote}" =~ /L([0-9]{3})/ ]]; then
    id="${BASH_REMATCH[1]}"
  else
    echo "WARN: Cannot parse sample ID from ${remote}, skipping."
    continue
  fi

  fname="$(basename "${remote}")"
  local_name="${fname}"
  if [[ "${id}" == "041" && "${remote}" == *"/P/"* ]]; then
    local_name="${fname/_P_/_N_}"
  fi
  local="${WORKDIR}/data/proteomics/l${id}/${local_name}"
  TASKS+=("${remote}"$'\t'"${local}")
done <<'EOF'
/mnt/data/lumanman/HCC_MS/L041/P/CNHPP_HCC_LC_profiling_L041_P_F1.raw
/mnt/data/lumanman/HCC_MS/L041/P/CNHPP_HCC_LC_profiling_L041_P_F2.raw
/mnt/data/lumanman/HCC_MS/L041/P/CNHPP_HCC_LC_profiling_L041_P_F3.raw
/mnt/data/lumanman/HCC_MS/L041/P/CNHPP_HCC_LC_profiling_L041_P_F4.raw
/mnt/data/lumanman/HCC_MS/L041/P/CNHPP_HCC_LC_profiling_L041_P_F5.raw
/mnt/data/lumanman/HCC_MS/L041/P/CNHPP_HCC_LC_profiling_L041_P_F6.raw
/mnt/data/lumanman/HCC_MS/L041/T/ms/CNHPP_HCC_LC_profiling_L041_T_F1.raw
/mnt/data/lumanman/HCC_MS/L041/T/ms/CNHPP_HCC_LC_profiling_L041_T_F2.raw
/mnt/data/lumanman/HCC_MS/L041/T/ms/CNHPP_HCC_LC_profiling_L041_T_F3.raw
/mnt/data/lumanman/HCC_MS/L041/T/ms/CNHPP_HCC_LC_profiling_L041_T_F4.raw
/mnt/data/lumanman/HCC_MS/L041/T/ms/CNHPP_HCC_LC_profiling_L041_T_F5.raw
/mnt/data/lumanman/HCC_MS/L041/T/ms/CNHPP_HCC_LC_profiling_L041_T_F6.raw
/mnt/data/lumanman/HCC_MS/L048/P/CNHPP_HCC_LC_profiling_L048P_5per_F1_R1.raw
/mnt/data/lumanman/HCC_MS/L048/P/CNHPP_HCC_LC_profiling_L048P_5per_F2_R1.raw
/mnt/data/lumanman/HCC_MS/L048/P/CNHPP_HCC_LC_profiling_L048P_5per_F3_R1.raw
/mnt/data/lumanman/HCC_MS/L048/P/CNHPP_HCC_LC_profiling_L048P_5per_F4_R1.raw
/mnt/data/lumanman/HCC_MS/L048/P/CNHPP_HCC_LC_profiling_L048P_5per_F5_R1.raw
/mnt/data/lumanman/HCC_MS/L048/P/CNHPP_HCC_LC_profiling_L048P_5per_F6_R1.raw
/mnt/data/lumanman/HCC_MS/L048/T/CNHPP_HCC_LC_profiling_L048T_5per_F1_R1.raw
/mnt/data/lumanman/HCC_MS/L048/T/CNHPP_HCC_LC_profiling_L048T_5per_F2_R1.raw
/mnt/data/lumanman/HCC_MS/L048/T/CNHPP_HCC_LC_profiling_L048T_5per_F3_R1.raw
/mnt/data/lumanman/HCC_MS/L048/T/CNHPP_HCC_LC_profiling_L048T_5per_F4_R1.raw
/mnt/data/lumanman/HCC_MS/L048/T/CNHPP_HCC_LC_profiling_L048T_5per_F5_R1.raw
/mnt/data/lumanman/HCC_MS/L048/T/CNHPP_HCC_LC_profiling_L048T_5per_F6_R1.raw
/mnt/data/lumanman/HCC_MS/L052/P/CNHPP_HCC_LC_profiling_L052P_5per_F1_R1.raw
/mnt/data/lumanman/HCC_MS/L052/P/CNHPP_HCC_LC_profiling_L052P_5per_F2_R1.raw
/mnt/data/lumanman/HCC_MS/L052/P/CNHPP_HCC_LC_profiling_L052P_5per_F3_R1.raw
/mnt/data/lumanman/HCC_MS/L052/P/CNHPP_HCC_LC_profiling_L052P_5per_F4_R1.raw
/mnt/data/lumanman/HCC_MS/L052/P/CNHPP_HCC_LC_profiling_L052P_5per_F5_R1.raw
/mnt/data/lumanman/HCC_MS/L052/P/CNHPP_HCC_LC_profiling_L052P_5per_F6_R1.raw
/mnt/data/lumanman/HCC_MS/L052/T/CNHPP_HCC_LC_profiling_L052T_5per_F1_R1.raw
/mnt/data/lumanman/HCC_MS/L052/T/CNHPP_HCC_LC_profiling_L052T_5per_F2_R1.raw
/mnt/data/lumanman/HCC_MS/L052/T/CNHPP_HCC_LC_profiling_L052T_5per_F3_R1.raw
/mnt/data/lumanman/HCC_MS/L052/T/CNHPP_HCC_LC_profiling_L052T_5per_F4_R1.raw
/mnt/data/lumanman/HCC_MS/L052/T/CNHPP_HCC_LC_profiling_L052T_5per_F5_R1.raw
/mnt/data/lumanman/HCC_MS/L052/T/CNHPP_HCC_LC_profiling_L052T_5per_F6_R1.raw
/mnt/data/lumanman/HCC_MS/L056/P/CNHPP_HCC_LC_profiling_L056P_5per_F1_R2.raw
/mnt/data/lumanman/HCC_MS/L056/P/CNHPP_HCC_LC_profiling_L056P_5per_F2_R2.raw
/mnt/data/lumanman/HCC_MS/L056/P/CNHPP_HCC_LC_profiling_L056P_5per_F3_R2.raw
/mnt/data/lumanman/HCC_MS/L056/P/CNHPP_HCC_LC_profiling_L056P_5per_F4_R2.raw
/mnt/data/lumanman/HCC_MS/L056/P/CNHPP_HCC_LC_profiling_L056P_5per_F5_R2.raw
/mnt/data/lumanman/HCC_MS/L056/P/CNHPP_HCC_LC_profiling_L056P_5per_F6_R1.raw
/mnt/data/lumanman/HCC_MS/L056/T/CNHPP_HCC_LC_profiling_L056T_5per_F1_R2.raw
/mnt/data/lumanman/HCC_MS/L056/T/CNHPP_HCC_LC_profiling_L056T_5per_F2_R2.raw
/mnt/data/lumanman/HCC_MS/L056/T/CNHPP_HCC_LC_profiling_L056T_5per_F3_R2.raw
/mnt/data/lumanman/HCC_MS/L056/T/CNHPP_HCC_LC_profiling_L056T_5per_F4_R2.raw
/mnt/data/lumanman/HCC_MS/L056/T/CNHPP_HCC_LC_profiling_L056T_5per_F5_R2.raw
/mnt/data/lumanman/HCC_MS/L056/T/CNHPP_HCC_LC_profiling_L056T_5per_F6_R1.raw
EOF

printf '%s\0' "${TASKS[@]}" | xargs -0 -P "${PARALLEL}" -I {} bash -c '
  line="$1"
  IFS=$'"'"'\t'"'"' read -r remote local <<< "$line"
  sync_one "$remote" "$local"
' _ {}
