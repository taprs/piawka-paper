#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 8 ]]; then
  echo "Usage: $0 <vcf.gz> <stats_csv> <mode:full|windowed> <window_or_NA> <out_prefix> <groups_tsv> <chrom_limit> <site_limit>" >&2
  exit 1
fi

VCF="$1"
STATS="$2"
MODE="$3"
WINDOW="$4"
OUT_PREFIX="$5"
POP_FILE="$6"
CHROM_LIMIT="$7"
SITE_LIMIT="$8"
PIXY_BIN="${PIXY_BIN:-pixy}"
PIXY_N_CORES="${PIXY_N_CORES:-}"

mkdir -p "$(dirname "${OUT_PREFIX}")"
WORKDIR="$(dirname "${OUT_PREFIX}")/pixy_tmp_$(basename "${OUT_PREFIX}")"
mkdir -p "${WORKDIR}"
TIME_LOG="${OUT_PREFIX}.time.tsv"
RAW_TIME_LOG="${OUT_PREFIX}.time.raw.log"

OUT_DIR="${OUT_PREFIX}.pixy"
mkdir -p "${OUT_DIR}"

# pixy population file must not contain a header row.
PIX_POP_FILE="${WORKDIR}/pixy_populations.tsv"
awk 'BEGIN{FS=OFS="\t"} NR==1 {h=tolower($1); if (h=="sample_name" || h=="sample" || h=="sampleid") next} NF>=2 {print $1,$2}' "${POP_FILE}" > "${PIX_POP_FILE}"

# Optional region downsampling for pilot runs using the first contig and a site cap.
INPUT_VCF="${VCF}"
if [[ "${CHROM_LIMIT}" != "NA" && "${SITE_LIMIT}" != "NA" ]]; then
  PILOT_VCF="${WORKDIR}/pilot_slice.vcf.gz"
  bcftools view -h "${VCF}" > "${WORKDIR}/header.vcf"
  bcftools view -H -r "${CHROM_LIMIT}" "${VCF}" | awk -v n="${SITE_LIMIT}" 'NR<=n' > "${WORKDIR}/body.vcf"
  cat "${WORKDIR}/header.vcf" "${WORKDIR}/body.vcf" | bgzip -c > "${PILOT_VCF}"
  tabix -f "${PILOT_VCF}"
  INPUT_VCF="${PILOT_VCF}"
fi

IFS=',' read -ra STATS_ARR <<< "${STATS}"
pixy_cmd=("${PIXY_BIN}" --vcf "${INPUT_VCF}" --populations "${PIX_POP_FILE}" --stats "${STATS_ARR[@]}" --output_folder "${OUT_DIR}" --output_prefix "$(basename "${OUT_PREFIX}")")
if [[ -n "${PIXY_N_CORES}" ]]; then
  pixy_cmd+=(--n_cores "${PIXY_N_CORES}")
fi
if [[ "${MODE}" == "windowed" && "${WINDOW}" != "NA" ]]; then
  if [[ -f "${WINDOW}" ]]; then
    pixy_cmd+=(--bed_file "${WINDOW}")
  else
    pixy_cmd+=(--window_size "${WINDOW}")
  fi
else
  # pixy requires window_size or bed_file; use a moderate full-run window to avoid NA-casting failures.
  pixy_cmd+=(--window_size 10000)
fi
pixy_cmd+=(--bypass_invariant_check)

if [[ "${DEBUG_PIXY_CMD:-0}" == "1" ]]; then
  printf 'POP_FILE=%s\n' "${POP_FILE}" >&2
  printf 'PIXY_CMD=' >&2
  printf '%q ' "${pixy_cmd[@]}" >&2
  printf '\n' >&2
fi

start_ns="$(date +%s%N)"
group_metrics_file="${RAW_TIME_LOG}.group_metrics.tsv"
if command -v setsid >/dev/null 2>&1; then
  setsid /usr/bin/time -l "${pixy_cmd[@]}" 2> "${RAW_TIME_LOG}" &
else
  /usr/bin/time -l "${pixy_cmd[@]}" 2> "${RAW_TIME_LOG}" &
fi
runner_pid=$!
(
  awk -v runner="${runner_pid}" -v out="${group_metrics_file}" '
    function t2s(t,    n,a,s,m,h) {
      gsub(/^ +| +$/, "", t)
      gsub(",", ".", t)
      n = split(t, a, ":")
      if (n == 2) { m = a[1] + 0; s = a[2] + 0; return m * 60 + s }
      if (n == 3) { h = a[1] + 0; m = a[2] + 0; s = a[3] + 0; return h * 3600 + m * 60 + s }
      return t + 0
    }
    function is_desc(pid,    parent) {
      if (pid == runner) return 1
      if (!(pid in ppid)) return 0
      if (pid in memo) return memo[pid]
      parent = ppid[pid]
      if (parent == pid || parent == 0 || parent == "") return memo[pid] = 0
      return memo[pid] = is_desc(parent)
    }
    function sample(    cmd,line,n,f,pid,rss,u,s,rss_sum,k) {
      delete ppid
      delete rssv
      delete uv
      delete sv
      delete memo
      cmd = "ps -axo pid=,ppid=,rss=,utime=,stime="
      rss_sum = 0
      while ((cmd | getline line) > 0) {
        gsub(/^ +/, "", line)
        n = split(line, f, /[[:space:]]+/)
        if (n < 5) continue
        pid = f[1]
        ppid[pid] = f[2]
        rssv[pid] = f[3] + 0
        uv[pid] = t2s(f[4])
        sv[pid] = t2s(f[5])
      }
      close(cmd)
      for (k in ppid) {
        if (!is_desc(k)) continue
        rss = rssv[k]
        u = uv[k]
        s = sv[k]
        rss_sum += rss
        if (u > u_max[k]) u_max[k] = u
        if (s > s_max[k]) s_max[k] = s
      }
      if (rss_sum > peak_rss) peak_rss = rss_sum
    }
    BEGIN {
      peak_rss = 0
      while (system("kill -0 " runner " >/dev/null 2>&1") == 0) {
        sample()
        system("sleep 0.2")
      }
      sample()
      total_u = 0
      total_s = 0
      for (pid in u_max) total_u += u_max[pid]
      for (pid in s_max) total_s += s_max[pid]
      printf "peak_group_rss_kb\t%.0f\n", peak_rss > out
      printf "group_cpu_user_sec\t%.3f\n", total_u >> out
      printf "group_cpu_sys_sec\t%.3f\n", total_s >> out
    }
  '
) &
monitor_pid=$!
set +e
wait "${runner_pid}"
cmd_status=$?
set -e
wait "${monitor_pid}" || true
if [[ ${cmd_status} -ne 0 ]]; then
  exit "${cmd_status}"
fi
end_ns="$(date +%s%N)"
elapsed_ms="$(( (end_ns - start_ns) / 1000000 ))"
elapsed_sec="$(LC_ALL=C awk -v ms="${elapsed_ms}" 'BEGIN { printf "%.3f", ms/1000 }')"
raw_rss="$(awk '/maximum resident set size/ { print $1; exit }' "${RAW_TIME_LOG}")"
if [[ -z "${raw_rss}" ]]; then
  main_rss_kb="NA"
elif [[ "${raw_rss}" =~ ^[0-9]+$ ]]; then
  # BSD time on macOS reports bytes; normalize to KB for downstream tables.
  main_rss_kb="$(( (raw_rss + 1023) / 1024 ))"
else
  main_rss_kb="NA"
fi
group_rss_kb="$(awk -F'\t' '$1=="peak_group_rss_kb"{print $2}' "${group_metrics_file}" 2>/dev/null || true)"
if [[ "${group_rss_kb}" =~ ^[0-9]+$ ]] && (( group_rss_kb > 0 )); then
  max_rss_kb="${group_rss_kb}"
elif [[ "${main_rss_kb}" =~ ^[0-9]+$ ]]; then
  max_rss_kb="${main_rss_kb}"
else
  max_rss_kb="NA"
fi

cpu_user_sec="$(awk '/real/ && /user/ && /sys/ {u=""; for(i=1;i<=NF;i++) if($i=="user") u=$(i-1); if(u!="") last=u} END{gsub(",",".",last); print last}' "${RAW_TIME_LOG}")"
cpu_sys_sec="$(awk '/real/ && /user/ && /sys/ {s=""; for(i=1;i<=NF;i++) if($i=="sys") s=$(i-1); if(s!="") last=s} END{gsub(",",".",last); print last}' "${RAW_TIME_LOG}")"
group_cpu_user_sec="$(awk -F'\t' '$1=="group_cpu_user_sec"{print $2}' "${group_metrics_file}" 2>/dev/null || true)"
group_cpu_sys_sec="$(awk -F'\t' '$1=="group_cpu_sys_sec"{print $2}' "${group_metrics_file}" 2>/dev/null || true)"
group_cpu_user_sec="${group_cpu_user_sec/,/.}"
group_cpu_sys_sec="${group_cpu_sys_sec/,/.}"
if [[ -n "${group_cpu_user_sec}" && "${group_cpu_user_sec}" != "NA" && "${cpu_user_sec}" != "NA" ]]; then
  cpu_user_sec="$(LC_ALL=C awk -v a="${cpu_user_sec}" -v b="${group_cpu_user_sec}" 'BEGIN{print (a+0>b+0)?a:b}')"
elif [[ -n "${group_cpu_user_sec}" && "${group_cpu_user_sec}" != "NA" ]]; then
  cpu_user_sec="${group_cpu_user_sec}"
fi
if [[ -n "${group_cpu_sys_sec}" && "${group_cpu_sys_sec}" != "NA" && "${cpu_sys_sec}" != "NA" ]]; then
  cpu_sys_sec="$(LC_ALL=C awk -v a="${cpu_sys_sec}" -v b="${group_cpu_sys_sec}" 'BEGIN{print (a+0>b+0)?a:b}')"
elif [[ -n "${group_cpu_sys_sec}" && "${group_cpu_sys_sec}" != "NA" ]]; then
  cpu_sys_sec="${group_cpu_sys_sec}"
fi
if [[ -z "${cpu_user_sec}" ]]; then
  cpu_user_sec="NA"
fi
if [[ -z "${cpu_sys_sec}" ]]; then
  cpu_sys_sec="NA"
fi
if [[ "${cpu_user_sec}" != "NA" && "${cpu_sys_sec}" != "NA" ]]; then
  cpu_total_sec="$(LC_ALL=C awk -v u="${cpu_user_sec}" -v s="${cpu_sys_sec}" 'BEGIN { printf "%.3f", u+s }')"
else
  cpu_total_sec="NA"
fi

printf "tool\telapsed_sec\tcpu_user_sec\tcpu_sys_sec\tcpu_total_sec\tmax_rss_kb\npixy\t%s\t%s\t%s\t%s\t%s\n" \
  "${elapsed_sec}" "${cpu_user_sec}" "${cpu_sys_sec}" "${cpu_total_sec}" "${max_rss_kb}" > "${TIME_LOG}"

echo "Wrote ${OUT_DIR}"
