#!/usr/bin/env bash
# Data gain and multiallelic pi (Figure 2, Supplementary Figure S2), real A. lyrata data plus the 16 simulated
# neutral replicates, piawka vs pixy:
#   1. Data gain (real data): pi in 10 kbp windows with both tools                          -> Figure 2A, 2B
#   2. Multiallelic pi (real data): piawka --mult and pixy --include_multiallelic_snps      -> Figure 2C, S2
#   3. Multiallelic pi (simulated, in-text claim): piawka --mult on the 16 neutral replicates, compared with the
#      biallelic runs of scripts/run_resource_sim16.sh (results/resource_sim16/<Rep>/win10000000_p1) -> run that first.
# Resource usage (Figure 1, Supplementary Figure S1) is scripts/run_resource_sim16.sh and scripts/run_a2_real_scaling.sh.
# Figure 3 (North/South gene-wise Fst) is separate; see scripts/run_northsouth_fst_comparison.sh.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_DIR="${ROOT_DIR}/results"
STAGE_DIR="${BENCH_DIR}/staging"
STAGE_SIM_DIR="${STAGE_DIR}/sim16"
STAGE_REAL_DIR="${STAGE_DIR}/real"
# RESULTS_SUFFIX (optional) sends every output to results/<suffix>/ instead of results/, so a rerun with other
# tool versions does not overwrite the existing runs; staging stays shared. analyze_benchmarks.py honours it too.
OUT_BASE="${BENCH_DIR}${RESULTS_SUFFIX:+/${RESULTS_SUFFIX}}"
GAIN_DIR="${OUT_BASE}/data_gain"
MULTI_DIR="${OUT_BASE}/multiallelic_check"
LOG_DIR="${OUT_BASE}/logs"
mkdir -p "${STAGE_REAL_DIR}" "${GAIN_DIR}" "${MULTI_DIR}" "${LOG_DIR}" "${OUT_BASE}/tables"

# Staged data files.
ONEPOP_LIST="${STAGE_SIM_DIR}/onepop.txt"                 # paths of the 16 neutral one-population VCFs
SYN_WIN="${STAGE_SIM_DIR}/onepop_win.bed"                 # one 10 Mbp window (whole file)
SYN_GROUPS="${STAGE_SIM_DIR}/onepop_groups.tsv"
LYRATA_VCF="${STAGE_REAL_DIR}/All_lyrata_final_allpos_1mbp.vcf.gz"
LYRATA_WIN="${STAGE_REAL_DIR}/lyrata_win.bed"
LYRATA_GROUPS="${STAGE_REAL_DIR}/lyrata_groups.tsv"

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
for f in "${ONEPOP_LIST}" "${SYN_WIN}" "${SYN_GROUPS}" "${LYRATA_WIN}" "${LYRATA_GROUPS}" "${LYRATA_VCF}"; do
  if [[ ! -f "${f}" ]]; then
    echo "Missing required input: ${f}" >&2
    exit 1
  fi
done
awk 'BEGIN{FS=OFS="\t"} NR==1 {h=tolower($1); if (h=="sample_name" || h=="sample" || h=="sampleid") next} NF>=2 {print $1,$2}' "${LYRATA_GROUPS}" > "${LOCAL_LYRATA_GROUPS_CLEAN}"

# --- Benchmark 1: data gain on real data (Figure 2A, 2B) ---
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
  if [[ ! -f "${REAL_DIR}/piawka.time.tsv" ]]; then
    "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${LYRATA_VCF}" "pi,lines" "windowed" "${LYRATA_WIN}" "${REAL_DIR}/piawka" "${LOCAL_LYRATA_GROUPS_CLEAN}"
  fi
  if [[ ! -f "${REAL_DIR}/pixy.time.tsv" ]]; then
    "${ROOT_DIR}/scripts/run_pixy_benchmark.sh" "${LYRATA_VCF}" "pi" "windowed" "${LYRATA_WIN}" "${REAL_DIR}/pixy" "${LOCAL_LYRATA_GROUPS_CLEAN}" "NA" "NA"
  fi
  printf "%s\n" "${current_hash}" > "${hash_stamp}"
}

# --- Benchmark 2/3: biallelic vs multiallelic pi (Figure 2C, S2; simulated in-text claim) ---
run_multiallelic_benchmark() {
  echo "Running multiallelic pi benchmark (piawka --mult)..."
  SIM_MULTI_DIR="${MULTI_DIR}/simulated_mult"
  REAL_MULTI_DIR="${MULTI_DIR}/real_mult"
  mkdir -p "${SIM_MULTI_DIR}" "${REAL_MULTI_DIR}"

  # Simulated: one --mult run per replicate; the biallelic piawka and pixy runs it is compared with are the
  # win10000000_p1 runs of scripts/run_resource_sim16.sh (analyze_benchmarks.py reads them from there).
  sim_jobs="${SIM_MULTI_DIR}/jobs.tsv"
  while IFS= read -r vcf; do
    [[ -z "${vcf}" || "${vcf}" =~ ^# ]] && continue
    [[ "${vcf}" = /* ]] || vcf="${ROOT_DIR}/${vcf}"   # list entries are relative to the project root
    printf '%s\t%s\n' "$(basename "${vcf}" .vcf.gz)" "${vcf}"
  done < "${ONEPOP_LIST}" > "${sim_jobs}"
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
