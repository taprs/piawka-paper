#!/usr/bin/env bash
# A2 (review response): scaling on the real A. lyrata VCF (1 Mbp slice, 1018 samples; 910 in the group table).
#  (1) sample-size series: nested random subsets of the 910 grouped samples (25..910), 1 and 8 processes, pi only
#  (2) group-count series: all 910 samples merged into 1, 2, 4, 8 or all 12 populations (largest kept separate, the rest merged), 1 process, pi only
# Runs sequentially (timings). Subset VCFs/groups are built outside the timed step.
# Versions: piawka 27e6ae3 (0.9.1 + 2 commits) (tools/piawka-27e6ae3), pixy 2.2.3 (env piawka-paper).
# Usage: scripts/run_a2_real_scaling.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ENVBIN="${ENVBIN:-/Users/ntikhomirov/mambaforge/envs/piawka-paper/bin}"
export PATH="${ROOT}/tools/piawka-27e6ae3:${ENVBIN}:${PATH}"
export TMPDIR=/tmp
REAL="${ROOT}/results/staging/real"
VCF="${REAL}/All_lyrata_final_allpos_1mbp.vcf.gz"
WIN="${REAL}/lyrata_win.bed"
GROUPS_ALL="${REAL}/lyrata_groups.noheader.tsv"
OUT="${ROOT}/results/a2_real"
mkdir -p "${OUT}/inputs"

# nested subsets: one fixed shuffle (seed 20261031), first n samples
SHUF="${OUT}/inputs/samples_shuffled.txt"
if [[ ! -s "${SHUF}" ]]; then
  awk -F'\t' '{print $1}' "${GROUPS_ALL}" | sort | awk 'BEGIN{srand(20261031)} {print rand() "\t" $0}' | sort -k1,1n | awk -F'\t' '{print $2}' > "${SHUF}"
fi

run_pair() {  # tag vcf groups procs
  tag="$1"; vcf="$2"; grp="$3"; procs="$4"; d="${OUT}/${tag}_p${procs}"
  mkdir -p "${d}"
  [[ -s "${d}/done" ]] && return 0
  PIAWKA_JOBS="${procs}" "${ROOT}/scripts/run_piawka_benchmark.sh" "${vcf}" "pi" windowed "${WIN}" "${d}/piawka" "${grp}"
  PIXY_N_CORES="${procs}" "${ROOT}/scripts/run_pixy_benchmark.sh" "${vcf}" "pi" windowed "${WIN}" "${d}/pixy" "${grp}" NA NA
  date > "${d}/done"
}

for n in 25 50 100 200 400 910; do
  in="${OUT}/inputs/n${n}"
  if [[ ! -s "${in}.vcf.gz.tbi" ]]; then
    head -n "${n}" "${SHUF}" > "${in}.samples"
    awk -F'\t' 'NR==FNR{k[$1]=1; next} ($1 in k)' "${in}.samples" "${GROUPS_ALL}" > "${in}.groups.tsv"
    bcftools view -S "${in}.samples" -Oz -o "${in}.vcf.gz" "${VCF}"; tabix -f -p vcf "${in}.vcf.gz"
  fi
  for p in 1 8; do run_pair "n${n}" "${in}.vcf.gz" "${in}.groups.tsv" "${p}"; done
done

# group-count series at n=910: merge populations, largest first
IN="${OUT}/inputs/n910"
mk() {  # k -> groups file with k populations (largest populations kept separate, the rest merged)
  k="$1"
  awk -F'\t' '{print $2}' "${GROUPS_ALL}" | sort | uniq -c | sort -k1,1nr | sed 's/^ *[0-9]* //' > "${OUT}/inputs/group_order.txt"
  awk -F'\t' -v k="${k}" 'NR==FNR{r[$0]=NR; next} {i=r[$2]; print $1 "\tG" (i<k?i:k)}' "${OUT}/inputs/group_order.txt" "${GROUPS_ALL}"
}
for k in 1 2 4 8 12; do
  g="${OUT}/inputs/groups_k${k}.tsv"
  [[ -s "${g}" ]] || mk "${k}" > "${g}"
  run_pair "k${k}" "${IN}.vcf.gz" "${g}" 1
done
echo "A2 real scaling finished"
