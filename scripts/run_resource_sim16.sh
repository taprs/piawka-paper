#!/usr/bin/env bash
# Resource benchmark on the 16 rescaled neutral one-population replicates (results/staging/sim16/onepop), pi only.
# Sequential (timings/memory), skips runs that already have a `done` stamp. Output: results/resource_sim16/<Rep>/<tag>/
#  (1) window-size sweep, 1 process: 10 Mbp cut into windows of 0.625, 1.25, 2.5, 5, 10 Mbp          -> Figure 1A
#  (2) process series, 1/2/4/8 processes on eight 1.25 Mbp windows (dispatched to processes)          -> Figure 1B
#      (the 1 process point is win1250000_p1 of the sweep, so tags win1250000_p{1,2,4,8})
#  (3) sample-size series, 10/25/50/75 samples (nested random subsets, monomorphic sites dropped), 1 and 8 processes
#      on eight 1.25 Mbp windows; 100 samples = win1250000_p{1,8}                                       -> Supp. Figure S1
# Versions: piawka 27e6ae3 (0.9.1 + 2 commits) (tools/piawka-27e6ae3), pixy 2.2.3 (env piawka-paper).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ENVBIN="${ENVBIN:-/Users/ntikhomirov/mambaforge/envs/piawka-paper/bin}"
export PATH="${ROOT}/tools/piawka-27e6ae3:${ENVBIN}:${PATH}"
export TMPDIR=/tmp
SIM="${ROOT}/results/staging/sim16"
OUT="${ROOT}/results/resource_sim16"
GROUPS_ALL="${SIM}/onepop_groups.tsv"
SAMPLE_N="${SAMPLE_N:-10 25 50 75}"
PROCS="${PROCS:-1 2 4 8}"
WINDOWS="${WINDOWS:-625000 1250000 2500000 5000000 10000000}"
mkdir -p "${OUT}/inputs"

for w in ${WINDOWS}; do
  awk -v w="${w}" -v OFS='\t' 'BEGIN{for (s=0; s<10000000; s+=w) print 1, s, s+w}' > "${OUT}/inputs/win_${w}.bed"
done
WIN8="${OUT}/inputs/win_1250000.bed"

run_pair() {  # dir vcf groups procs bed
  local d="$1" vcf="$2" grp="$3" procs="$4" bed="$5"
  mkdir -p "${d}"; [[ -s "${d}/done" ]] && return 0
  PIAWKA_JOBS="${procs}" "${ROOT}/scripts/run_piawka_benchmark.sh" "${vcf}" pi windowed "${bed}" "${d}/piawka" "${grp}"
  PIXY_N_CORES="${procs}" "${ROOT}/scripts/run_pixy_benchmark.sh" "${vcf}" pi windowed "${bed}" "${d}/pixy" "${grp}" NA NA
  bcftools view -H "${vcf}" | wc -l > "${d}/records.txt"
  date > "${d}/done"
}

for r in $(seq 0 15); do
  rep="Demoonepop_scaled10_Rep${r}"
  vcf="${SIM}/onepop/${rep}.vcf.gz"
  [[ -f "${vcf}" ]] || { echo "missing ${vcf}" >&2; exit 1; }
  echo "== ${rep}"
  for w in ${WINDOWS}; do run_pair "${OUT}/${rep}/win${w}_p1" "${vcf}" "${GROUPS_ALL}" 1 "${OUT}/inputs/win_${w}.bed"; done
  for p in ${PROCS}; do [[ "${p}" == 1 ]] || run_pair "${OUT}/${rep}/win1250000_p${p}" "${vcf}" "${GROUPS_ALL}" "${p}" "${WIN8}"; done
  # nested subsets: first n of one shuffle of the sample names (seeded per replicate)
  shuf="${OUT}/inputs/${rep}.shuf"
  [[ -s "${shuf}" ]] || bcftools query -l "${vcf}" | awk -v s=$((20261108 + r)) 'BEGIN{srand(s)} {print rand() "\t" $0}' | sort -k1,1n | awk -F'\t' '{print $2}' > "${shuf}"
  for n in ${SAMPLE_N}; do
    in="${OUT}/inputs/${rep}_n${n}"
    if [[ ! -s "${in}.vcf.gz.tbi" ]]; then
      head -n "${n}" "${shuf}" > "${in}.samples"
      awk -v OFS='\t' '{print $1,"POP1"}' "${in}.samples" > "${in}.groups.tsv"
      bcftools view -S "${in}.samples" -c1 -Oz -o "${in}.vcf.gz" "${vcf}"; tabix -f -p vcf "${in}.vcf.gz"
    fi
    for p in 1 8; do run_pair "${OUT}/${rep}/n${n}_p${p}" "${in}.vcf.gz" "${in}.groups.tsv" "${p}" "${WIN8}"; done
    rm -f "${in}.vcf.gz" "${in}.vcf.gz.tbi"   # subsets are cheap to regenerate; keep the disk free
  done
done
echo "resource_sim16 finished"
