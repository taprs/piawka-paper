#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_DIR="${ROOT_DIR}/results"
INPUT_DIR="${BENCH_DIR}/inputs"
STAGE_DIR="${BENCH_DIR}/staging"
STAGE_SIM_DIR="${STAGE_DIR}/simulated"
STAGE_REAL_DIR="${STAGE_DIR}/real"
RES_DIR="${BENCH_DIR}/resource_usage"
ACC_DIR="${BENCH_DIR}/accuracy"
GAIN_DIR="${BENCH_DIR}/data_gain"
MULTI_DIR="${BENCH_DIR}/multiallelic_check"
LOG_DIR="${BENCH_DIR}/logs"
mkdir -p "${INPUT_DIR}" "${STAGE_SIM_DIR}" "${STAGE_REAL_DIR}" "${RES_DIR}" "${ACC_DIR}" "${GAIN_DIR}" "${MULTI_DIR}" "${LOG_DIR}" "${BENCH_DIR}/tables"

# Updated instruction paths.
NIKITA_ROOT="/Users/ntikhomirov/MPIPZ/netscratch/dep_mercier/grp_novikova/nikita/piawka_paper"
ONEPOP_LIST="${NIKITA_ROOT}/onepop.txt"
BNECK_LIST="${NIKITA_ROOT}/bneck.txt"
SYN_WIN="${NIKITA_ROOT}/syn_win.bed"
SYN_GROUPS="${NIKITA_ROOT}/syn_groups.tsv"
LYRATA_VCF="${NIKITA_ROOT}/All_lyrata_final_allpos_1mbp.vcf.gz"
LYRATA_WIN="${NIKITA_ROOT}/lyrata_win.bed"
LYRATA_GROUPS="${NIKITA_ROOT}/lyrata_groups.tsv"
ADD_MISSING_AWK="/Users/ntikhomirov/MPIPZ/netscratch/dep_mercier/grp_novikova/nikita/scripts/vcf_analysis/add_missing.awk"
MISSING_RATES=("0.05" "0.10" "0.20" "0.40")
RESOURCE_THREADS=("1" "2" "4" "8")

file_hash() {
  local f="$1"
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${f}" | awk '{print $1}'
  else
    md5 -q "${f}"
  fi
}

if [[ -x "${ROOT_DIR}/tools/piawka/piawka" ]]; then
  export PIAWKA_BIN="${ROOT_DIR}/tools/piawka/piawka"
fi
if [[ -x "/Users/ntikhomirov/mambaforge/envs/piawka-paper-py/bin/pixy" ]]; then
  export PIXY_BIN="/Users/ntikhomirov/mambaforge/envs/piawka-paper-py/bin/pixy"
fi

# Copy all static benchmark inputs into workspace to reduce filesystem latency.
LOCAL_SYN_WIN="${STAGE_SIM_DIR}/syn_win.bed"
LOCAL_SYN_GROUPS="${STAGE_SIM_DIR}/syn_groups.tsv"
LOCAL_LYRATA_WIN="${STAGE_REAL_DIR}/lyrata_win.bed"
LOCAL_LYRATA_GROUPS="${STAGE_REAL_DIR}/lyrata_groups.tsv"
LOCAL_LYRATA_GROUPS_CLEAN="${STAGE_REAL_DIR}/lyrata_groups.noheader.tsv"
LOCAL_LYRATA_VCF="${STAGE_REAL_DIR}/$(basename "${LYRATA_VCF}")"
LOCAL_LYRATA_TBI="${LOCAL_LYRATA_VCF}.tbi"
ensure_local_file() {
  local src="$1"
  local dst="$2"
  local label="$3"
  if [[ -f "${src}" ]]; then
    cp -f "${src}" "${dst}"
  elif [[ ! -f "${dst}" ]]; then
    echo "Missing required input: ${label} (${src} or ${dst})" >&2
    exit 1
  fi
}
ensure_local_file "${SYN_WIN}" "${LOCAL_SYN_WIN}" "syn window bed"
ensure_local_file "${SYN_GROUPS}" "${LOCAL_SYN_GROUPS}" "syn groups file"
ensure_local_file "${LYRATA_WIN}" "${LOCAL_LYRATA_WIN}" "lyrata window bed"
ensure_local_file "${LYRATA_GROUPS}" "${LOCAL_LYRATA_GROUPS}" "lyrata groups file"
ensure_local_file "${LYRATA_VCF}" "${LOCAL_LYRATA_VCF}" "lyrata vcf"
awk 'BEGIN{FS=OFS="\t"} NR==1 {h=tolower($1); if (h=="sample_name" || h=="sample" || h=="sampleid") next} NF>=2 {print $1,$2}' "${LOCAL_LYRATA_GROUPS}" > "${LOCAL_LYRATA_GROUPS_CLEAN}"
if [[ -f "${LYRATA_VCF}.tbi" && ! -f "${LOCAL_LYRATA_TBI}" ]]; then
  cp -f "${LYRATA_VCF}.tbi" "${LOCAL_LYRATA_TBI}"
fi

# Use a workspace-local copy of add_missing.awk to avoid intermittent external FS failures.
LOCAL_ADD_MISSING_AWK="${INPUT_DIR}/add_missing.awk"
if [[ -f "${ADD_MISSING_AWK}" ]]; then
  cp "${ADD_MISSING_AWK}" "${LOCAL_ADD_MISSING_AWK}"
elif [[ ! -f "${LOCAL_ADD_MISSING_AWK}" ]]; then
  echo "Missing required input: add_missing.awk (${ADD_MISSING_AWK} or ${LOCAL_ADD_MISSING_AWK})" >&2
  exit 1
fi

MANIFEST_NOMISS="${INPUT_DIR}/sim_nomissing_manifest.tsv"
MANIFEST_MISSING="${INPUT_DIR}/sim_missing_manifest.tsv"
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

echo -e "sample_id\tsource_class\tmissing_rate\tvcf_path" > "${MANIFEST_MISSING}"

echo "Preparing missing-data simulated inputs..."
tail -n +2 "${MANIFEST_NOMISS}" | while IFS=$'\t' read -r sample_id source_class vcf_path; do
  for rate in "${MISSING_RATES[@]}"; do
    rate_tag="r$(echo "${rate}" | tr -d '.')"
    out_dir="${INPUT_DIR}/missing/${rate_tag}"
    out_vcf="${out_dir}/${sample_id}.vcf.gz"
    mkdir -p "${out_dir}"
    if [[ ! -f "${out_vcf}" || ! -f "${out_vcf}.tbi" ]]; then
      bgzip -cd "${vcf_path}" | awk -v MIN="${rate}" -v MAX="${rate}" -f "${LOCAL_ADD_MISSING_AWK}" | bgzip -c > "${out_vcf}"
      tabix -f "${out_vcf}"
    fi
    echo -e "${sample_id}\t${source_class}\t${rate}\t${out_vcf}" >> "${MANIFEST_MISSING}"
  done
done

echo "Running resource usage benchmark (simulated, no missing)..."
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
      PIAWKA_JOBS="${threads}" "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${vcf_path}" "pi,lines,miss,theta_w" "windowed" "${LOCAL_SYN_WIN}" "${piawka_prefix}" "${LOCAL_SYN_GROUPS}"
    fi
    if [[ ! -f "${pixy_prefix}.time.tsv" ]]; then
      PIXY_N_CORES="${threads}" "${ROOT_DIR}/scripts/run_pixy_benchmark.sh" "${vcf_path}" "pi,watterson_theta" "windowed" "${LOCAL_SYN_WIN}" "${pixy_prefix}" "${LOCAL_SYN_GROUPS}" "NA" "NA"
    fi
  done
done

run_accuracy_benchmark() {
  echo "Running accuracy benchmark (simulated + missing data)..."
  tail -n +2 "${MANIFEST_MISSING}" | while IFS=$'\t' read -r sample_id _ rate vcf_path; do
    rate_tag="r$(echo "${rate}" | tr -d '.')"
    run_dir="${ACC_DIR}/${sample_id}/${rate_tag}"
    mkdir -p "${run_dir}"
    if [[ ! -f "${run_dir}/piawka.time.tsv" ]]; then
      "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${vcf_path}" "pi,theta_w,theta_low" "windowed" "${LOCAL_SYN_WIN}" "${run_dir}/piawka" "${LOCAL_SYN_GROUPS}"
    fi
    if [[ ! -f "${run_dir}/pixy.time.tsv" ]]; then
      "${ROOT_DIR}/scripts/run_pixy_benchmark.sh" "${vcf_path}" "pi,watterson_theta" "windowed" "${LOCAL_SYN_WIN}" "${run_dir}/pixy" "${LOCAL_SYN_GROUPS}" "NA" "NA"
    fi
  done
}

run_biallelic_check() {
  echo "Running biallelic-vs-multiallelic pi check..."
  first_sample="$(awk 'BEGIN{FS="\t"} NR==2 {print $1; exit}' "${MANIFEST_NOMISS}")"
  first_vcf="$(awk 'BEGIN{FS="\t"} NR==2 {print $3; exit}' "${MANIFEST_NOMISS}")"
  if [[ -n "${first_sample}" && -n "${first_vcf}" ]]; then
    bial_vcf="${MULTI_DIR}/${first_sample}.bial.vcf.gz"
    if [[ ! -f "${bial_vcf}" || ! -f "${bial_vcf}.tbi" ]]; then
      bcftools view -m2 -M2 -v snps "${first_vcf}" -Oz -o "${bial_vcf}"
      tabix -f "${bial_vcf}"
    fi
    "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${bial_vcf}" "pi" "windowed" "${LOCAL_SYN_WIN}" "${MULTI_DIR}/${first_sample}.piawka_bial" "${LOCAL_SYN_GROUPS}"
  fi
}

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
  "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${LOCAL_LYRATA_VCF}" "pi,lines" "windowed" "${LOCAL_LYRATA_WIN}" "${REAL_DIR}/piawka" "${LOCAL_LYRATA_GROUPS_CLEAN}"
  "${ROOT_DIR}/scripts/run_pixy_benchmark.sh" "${LOCAL_LYRATA_VCF}" "pi" "windowed" "${LOCAL_LYRATA_WIN}" "${REAL_DIR}/pixy" "${LOCAL_LYRATA_GROUPS_CLEAN}" "NA" "NA"
  printf "%s\n" "${current_hash}" > "${hash_stamp}"
}

run_multiallelic_benchmark() {
  echo "Running multiallelic pi benchmark (piawka --mult)..."
  SIM_MULTI_DIR="${MULTI_DIR}/simulated_mult"
  REAL_MULTI_DIR="${MULTI_DIR}/real_mult"
  mkdir -p "${SIM_MULTI_DIR}" "${REAL_MULTI_DIR}"

  sim_jobs="${SIM_MULTI_DIR}/jobs.tsv"
  tail -n +2 "${MANIFEST_NOMISS}" | awk -F'\t' '{print $1"\t"$3}' > "${sim_jobs}"
  parallel -j 8 --colsep '\t' \
    "if [[ ! -f '${SIM_MULTI_DIR}/{1}/piawka_mult.time.tsv' ]] || ! grep -q \$'\\tlines\\t' '${SIM_MULTI_DIR}/{1}/piawka_mult.piawka.bed' 2>/dev/null; then '${ROOT_DIR}/scripts/run_piawka_benchmark.sh' '{2}' 'pi,lines' 'windowed' '${LOCAL_SYN_WIN}' '${SIM_MULTI_DIR}/{1}/piawka_mult' '${LOCAL_SYN_GROUPS}' '--mult'; fi" \
    :::: "${sim_jobs}"

  real_hash_stamp="${REAL_MULTI_DIR}/groups.sha256"
  current_hash="$(file_hash "${LOCAL_LYRATA_GROUPS_CLEAN}")"
  previous_hash="$(cat "${real_hash_stamp}" 2>/dev/null || true)"
  if [[ "${current_hash}" != "${previous_hash}" ]]; then
    rm -f "${REAL_MULTI_DIR}/piawka_mult.piawka.bed" "${REAL_MULTI_DIR}/piawka_mult.time.tsv" "${REAL_MULTI_DIR}/piawka_mult.time.raw.log" "${REAL_MULTI_DIR}/piawka_mult.time.raw.log.group_metrics.tsv"
  fi
  if [[ ! -f "${REAL_MULTI_DIR}/piawka_mult.time.tsv" ]] || ! grep -q $'\tlines\t' "${REAL_MULTI_DIR}/piawka_mult.piawka.bed" 2>/dev/null; then
    "${ROOT_DIR}/scripts/run_piawka_benchmark.sh" "${LOCAL_LYRATA_VCF}" "pi,lines" "windowed" "${LOCAL_LYRATA_WIN}" "${REAL_MULTI_DIR}/piawka_mult" "${LOCAL_LYRATA_GROUPS_CLEAN}" "--mult"
  fi
  printf "%s\n" "${current_hash}" > "${real_hash_stamp}"
}

# Per request: run all benchmarks except resource usage in parallel.
run_accuracy_benchmark > "${LOG_DIR}/accuracy.log" 2>&1 &
pid_accuracy=$!
run_data_gain_benchmark > "${LOG_DIR}/data_gain.log" 2>&1 &
pid_data_gain=$!
run_biallelic_check > "${LOG_DIR}/biallelic_check.log" 2>&1 &
pid_bial=$!
run_multiallelic_benchmark > "${LOG_DIR}/multiallelic.log" 2>&1 &
pid_multi=$!
wait_status=0
for pid in "${pid_accuracy}" "${pid_data_gain}" "${pid_bial}" "${pid_multi}"; do
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
