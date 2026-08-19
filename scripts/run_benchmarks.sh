#!/usr/bin/env bash
# Runs the three benchmarks reported in the paper:
#   1. Resource usage (simulated data, piawka vs pixy at 1/2/4/8 parallel processes)  -> Figure 1B
#   2. Window-size effect on memory (simulated data, 1 process, 10 Mbp vs 8x1.25 Mbp) -> Figure 1A
#   3. Data gain + multiallelic pi (real A. lyrata data, and simulated for the
#      multiallelic comparison)                                                       -> Figure 2
#
# Figure 3 (North/South gene-wise Fst) is a separate, independent analysis; see
# scripts/run_northsouth_fst_comparison.sh.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_DIR="${ROOT_DIR}/results"
INPUT_DIR="${BENCH_DIR}/inputs"
STAGE_DIR="${BENCH_DIR}/staging"
STAGE_SIM_DIR="${STAGE_DIR}/simulated"
STAGE_REAL_DIR="${STAGE_DIR}/real"
RES_DIR="${BENCH_DIR}/resource_usage"
WIN_DIR="${BENCH_DIR}/resource_usage_compare_1thread"
GAIN_DIR="${BENCH_DIR}/data_gain"
MULTI_DIR="${BENCH_DIR}/multiallelic_check"
LOG_DIR="${BENCH_DIR}/logs"
mkdir -p "${INPUT_DIR}" "${STAGE_SIM_DIR}" "${STAGE_REAL_DIR}" "${RES_DIR}" "${WIN_DIR}" "${GAIN_DIR}" "${MULTI_DIR}" "${LOG_DIR}" "${BENCH_DIR}/tables"

# Staged data files (relative to project root).
ONEPOP_LIST="${STAGE_SIM_DIR}/onepop.txt"
BNECK_LIST="${STAGE_SIM_DIR}/bneck.txt"
SYN_WIN="${STAGE_SIM_DIR}/syn_win.bed"                    # one 10 Mbp window (whole file)
SYN_WIN_8="${STAGE_SIM_DIR}/syn_win_threads8.bed"         # eight equal 1.25 Mbp windows
SYN_GROUPS="${STAGE_SIM_DIR}/syn_groups.tsv"
LYRATA_VCF="${STAGE_REAL_DIR}/All_lyrata_final_allpos_1mbp.vcf.gz"
LYRATA_WIN="${STAGE_REAL_DIR}/lyrata_win.bed"
LYRATA_GROUPS="${STAGE_REAL_DIR}/lyrata_groups.tsv"
RESOURCE_THREADS=("1" "2" "4" "8")

file_hash() {
  local f="$1"
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${f}" | awk '{print $1}'
  else
    md5 -q "${f}"
  fi
}

if command -v pixy >/dev/null 2>&1; then
  export PIXY_BIN="$(command -v pixy)"
fi

LOCAL_LYRATA_GROUPS_CLEAN="${STAGE_REAL_DIR}/lyrata_groups.noheader.tsv"
for f in "${SYN_WIN}" "${SYN_GROUPS}" "${LYRATA_WIN}" "${LYRATA_GROUPS}" "${LYRATA_VCF}"; do
  if [[ ! -f "${f}" ]]; then
    echo "Missing required input: ${f}" >&2
    exit 1
  fi
done
awk 'BEGIN{FS=OFS="\t"} NR==1 {h=tolower($1); if (h=="sample_name" || h=="sample" || h=="sampleid") next} NF>=2 {print $1,$2}' "${LYRATA_GROUPS}" > "${LOCAL_LYRATA_GROUPS_CLEAN}"

# The 8-window BED is a plain split of the single-window BED; derive it if absent.
if [[ ! -f "${SYN_WIN_8}" ]]; then
  awk 'BEGIN{FS=OFS="\t"} {n=8; w=int(($3-$2)/n); for (i=0;i<n;i++) print $1, $2+i*w, (i==n-1 ? $3 : $2+(i+1)*w)}' "${SYN_WIN}" > "${SYN_WIN_8}"
fi

MANIFEST_NOMISS="${INPUT_DIR}/sim_nomissing_manifest.tsv"
emit_manifest_rows() {
  local class_name="$1"
  local list_file="$2"
  while IFS= read -r vcf; do
    [[ -z "${vcf}" ]] && continue
    [[ "${vcf}" =~ ^# ]] && continue
    # Normalize whitespace and leading-tilde paths from list files.
    vcf="$(echo "${vcf}" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    vcf="${vcf/#\~/${HOME}}"
    [[ -f "${vcf}" ]] || { echo "Skipping missing VCF: ${vcf}" >&2; continue; }
    local sample_id
    sample_id="$(basename "$(dirname "${vcf}")")"
    local local_vcf="${STAGE_SIM_DIR}/${sample_id}.vcf.gz"
    local local_tbi="${local_vcf}.tbi"
    if [[ ! -f "${local_vcf}" ]]; then
      cp -f "${vcf}" "${local_vcf}"
    fi
    if [[ -f "${vcf}.tbi" && ! -f "${local_tbi}" ]]; then
      cp -f "${vcf}.tbi" "${local_tbi}"
    fi
    echo -e "${sample_id}\t${class_name}\t${local_vcf}" >> "${MANIFEST_NOMISS}"
  done < "${list_file}"
}

if [[ -f "${ONEPOP_LIST}" && -f "${BNECK_LIST}" ]]; then
  echo -e "sample_id\tsource_class\tvcf_path" > "${MANIFEST_NOMISS}"
  emit_manifest_rows "onepop" "${ONEPOP_LIST}"
  emit_manifest_rows "bneck" "${BNECK_LIST}"
  # Remove duplicate sample IDs keeping first entry.
  awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} !seen[$1]++ {print}' "${MANIFEST_NOMISS}" > "${MANIFEST_NOMISS}.tmp" && mv "${MANIFEST_NOMISS}.tmp" "${MANIFEST_NOMISS}"
elif [[ -s "${MANIFEST_NOMISS}" ]]; then
  echo "Using existing simulated manifest: ${MANIFEST_NOMISS}"
else
  echo "Missing required inputs: onepop/bneck lists not found and no existing manifest (${MANIFEST_NOMISS})" >&2
  exit 1
fi

# --- Benchmark 1: resource usage vs number of parallel processes (Figure 1B) ---
echo "Running resource usage benchmark (simulated)..."
tail -n +2 "${MANIFEST_NOMISS}" | while IFS=$'\t' read -r sample_id _ vcf_path; do
  run_dir="${RES_DIR}/${sample_id}"
  mkdir -p "${run_dir}"
  for threads in "${RESOURCE_THREADS[@]}"; do
    suffix=""
    if [[ "${threads}" != "1" ]]; then
      suffix="_t${threads}"
    fi
    piawka_prefix="${run_dir}/piawka${suffix}"
    pixy_prefix="${run_dir}/pixy${suffix}"
    if [[ ! -f "${piawka_prefix}.time.tsv" ]]; then
      PIAWKA_JOBS="${threads}" "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${vcf_path}" "pi,lines,miss" "windowed" "${SYN_WIN}" "${piawka_prefix}" "${SYN_GROUPS}"
    fi
    if [[ ! -f "${pixy_prefix}.time.tsv" ]]; then
      PIXY_N_CORES="${threads}" "${ROOT_DIR}/scripts/run_pixy_benchmark.sh" "${vcf_path}" "pi" "windowed" "${SYN_WIN}" "${pixy_prefix}" "${SYN_GROUPS}" "NA" "NA"
    fi
  done
done

# --- Benchmark 2: window-size effect on memory, single process (Figure 1A) ---
# Same VCFs, same single process, two BED files: one 10 Mbp window vs eight
# 1.25 Mbp windows. Aggregated into
# results/tables/resource_usage_1thread_synwin_vs_threads8.tsv.
run_window_size_benchmark() {
  echo "Running window-size comparison benchmark (simulated, 1 process)..."
  tail -n +2 "${MANIFEST_NOMISS}" | while IFS=$'\t' read -r sample_id _ vcf_path; do
    run_dir="${WIN_DIR}/${sample_id}"
    mkdir -p "${run_dir}"
    for win_tag in syn_win syn_win_threads8; do
      if [[ "${win_tag}" == "syn_win" ]]; then
        bed="${SYN_WIN}"
      else
        bed="${SYN_WIN_8}"
      fi
      if [[ ! -f "${run_dir}/piawka_${win_tag}.time.tsv" ]]; then
        PIAWKA_JOBS="1" "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${vcf_path}" "pi,lines,miss" "windowed" "${bed}" "${run_dir}/piawka_${win_tag}" "${SYN_GROUPS}"
      fi
      if [[ ! -f "${run_dir}/pixy_${win_tag}.time.tsv" ]]; then
        PIXY_N_CORES="1" "${ROOT_DIR}/scripts/run_pixy_benchmark.sh" "${vcf_path}" "pi" "windowed" "${bed}" "${run_dir}/pixy_${win_tag}" "${SYN_GROUPS}" "NA" "NA"
      fi
    done
  done
}

# --- Benchmark 3a: data gain on real data (Figure 2A, 2B) ---
run_data_gain_benchmark() {
  echo "Running data gain benchmark (real data windows)..."
  REAL_DIR="${GAIN_DIR}/real_lyrata"
  mkdir -p "${REAL_DIR}"
  hash_stamp="${REAL_DIR}/groups.sha256"
  current_hash="$(file_hash "${LOCAL_LYRATA_GROUPS_CLEAN}")"
  previous_hash="$(cat "${hash_stamp}" 2>/dev/null || true)"
  if [[ "${current_hash}" != "${previous_hash}" ]]; then
    rm -f "${REAL_DIR}/piawka.piawka.bed" "${REAL_DIR}/piawka.time.tsv" "${REAL_DIR}/piawka.time.raw.log" "${REAL_DIR}/piawka.time.raw.log.group_metrics.tsv"
    rm -f "${REAL_DIR}/pixy.time.tsv" "${REAL_DIR}/pixy.time.raw.log" "${REAL_DIR}/pixy.time.raw.log.group_metrics.tsv"
    rm -rf "${REAL_DIR}/pixy.pixy"
  fi
  "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${LYRATA_VCF}" "pi,lines" "windowed" "${LYRATA_WIN}" "${REAL_DIR}/piawka" "${LOCAL_LYRATA_GROUPS_CLEAN}"
  "${ROOT_DIR}/scripts/run_pixy_benchmark.sh" "${LYRATA_VCF}" "pi" "windowed" "${LYRATA_WIN}" "${REAL_DIR}/pixy" "${LOCAL_LYRATA_GROUPS_CLEAN}" "NA" "NA"
  printf "%s\n" "${current_hash}" > "${hash_stamp}"
}

# --- Benchmark 3b: biallelic vs multiallelic pi (Figure 2C, plus simulated) ---
run_multiallelic_benchmark() {
  echo "Running multiallelic pi benchmark (piawka --mult)..."
  SIM_MULTI_DIR="${MULTI_DIR}/simulated_mult"
  REAL_MULTI_DIR="${MULTI_DIR}/real_mult"
  mkdir -p "${SIM_MULTI_DIR}" "${REAL_MULTI_DIR}"

  sim_jobs="${SIM_MULTI_DIR}/jobs.tsv"
  tail -n +2 "${MANIFEST_NOMISS}" | awk -F'\t' '{print $1"\t"$3}' > "${sim_jobs}"
  parallel -j 8 --colsep '\t' \
    "if [[ ! -f '${SIM_MULTI_DIR}/{1}/piawka_mult.time.tsv' ]] || ! grep -q \$'\\tlines\\t' '${SIM_MULTI_DIR}/{1}/piawka_mult.piawka.bed' 2>/dev/null; then '${ROOT_DIR}/scripts/run_piawka_benchmark.sh' '{2}' 'pi,lines' 'windowed' '${SYN_WIN}' '${SIM_MULTI_DIR}/{1}/piawka_mult' '${SYN_GROUPS}' '--mult'; fi" \
    :::: "${sim_jobs}"

  real_hash_stamp="${REAL_MULTI_DIR}/groups.sha256"
  current_hash="$(file_hash "${LOCAL_LYRATA_GROUPS_CLEAN}")"
  previous_hash="$(cat "${real_hash_stamp}" 2>/dev/null || true)"
  if [[ "${current_hash}" != "${previous_hash}" ]]; then
    rm -f "${REAL_MULTI_DIR}/piawka_mult.piawka.bed" "${REAL_MULTI_DIR}/piawka_mult.time.tsv" "${REAL_MULTI_DIR}/piawka_mult.time.raw.log" "${REAL_MULTI_DIR}/piawka_mult.time.raw.log.group_metrics.tsv"
  fi
  if [[ ! -f "${REAL_MULTI_DIR}/piawka_mult.time.tsv" ]] || ! grep -q $'\tlines\t' "${REAL_MULTI_DIR}/piawka_mult.piawka.bed" 2>/dev/null; then
    "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${LYRATA_VCF}" "pi,lines" "windowed" "${LYRATA_WIN}" "${REAL_MULTI_DIR}/piawka_mult" "${LOCAL_LYRATA_GROUPS_CLEAN}" "--mult"
  fi
  if [[ "${current_hash}" != "${previous_hash}" ]]; then
    rm -f "${REAL_MULTI_DIR}/pixy_mult.time.tsv" "${REAL_MULTI_DIR}/pixy_mult.time.raw.log" "${REAL_MULTI_DIR}/pixy_mult.time.raw.log.group_metrics.tsv"
    rm -rf "${REAL_MULTI_DIR}/pixy_mult.pixy"
  fi
  if [[ ! -f "${REAL_MULTI_DIR}/pixy_mult.time.tsv" ]]; then
    "${ROOT_DIR}/scripts/run_pixy_benchmark.sh" "${LYRATA_VCF}" "pi" "windowed" "${LYRATA_WIN}" "${REAL_MULTI_DIR}/pixy_mult" "${LOCAL_LYRATA_GROUPS_CLEAN}" "NA" "NA" "--include_multiallelic_snps"
  fi
  printf "%s\n" "${current_hash}" > "${real_hash_stamp}"
}

# The window-size benchmark measures memory and must not share the machine with
# the other stages; run it on its own, then the rest in parallel.
run_window_size_benchmark > "${LOG_DIR}/window_size.log" 2>&1

run_data_gain_benchmark > "${LOG_DIR}/data_gain.log" 2>&1 &
pid_data_gain=$!
run_multiallelic_benchmark > "${LOG_DIR}/multiallelic.log" 2>&1 &
pid_multi=$!
wait_status=0
for pid in "${pid_data_gain}" "${pid_multi}"; do
  if ! wait "${pid}"; then
    wait_status=1
  fi
done
if [[ "${wait_status}" -ne 0 ]]; then
  exit 1
fi

echo "Aggregating benchmark results..."
python3 "${ROOT_DIR}/scripts/analyze_benchmarks.py"
echo "Benchmark run complete."
