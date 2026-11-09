#!/usr/bin/env bash
# A2.4 (review response): where does piawka spend its time? Stage ablation on a 200 kb region of the real
# A. lyrata VCF (910 samples, 1 process, pi only unless stated). Each stage adds one layer; run alone on an idle machine.
#   a  tabix region -> /dev/null                         (decompression + I/O floor)
#   b  + gawk field splitting only                       (NF per record)
#   c  + genotype scan of every sample column            (the per-sample allele-counting loop of piawka's process_sites, no statistics)
#   d  piawka calc -s pi   (1 group)   e  piawka calc -s pi,theta_w,tajima,lines,miss   f  piawka calc -s pi, 12 groups (all population labels)
#   g  d under LC_ALL=en_US.UTF-8 (piawka itself exports LC_ALL=C)
# Also writes a gawk --profile execution-count listing for piawka on a 20 kb region.
# Usage: scripts/run_a2_profile.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ENVBIN="${ENVBIN:-/Users/ntikhomirov/mambaforge/envs/piawka-paper/bin}"
export PATH="${ROOT}/tools/piawka-27e6ae3:${ENVBIN}:${PATH}"
export TMPDIR=/tmp
IN="${ROOT}/results/a2_real/inputs"
VCF="${IN}/n910.vcf.gz"
OUT="${ROOT}/results/a2_profile"; mkdir -p "${OUT}"
REG="scaffold_1:1000001-1200000"
awk -v OFS='\t' '{print $1, $2, $3}' "${ROOT}/results/staging/real/lyrata_win.bed" | head -n 20 > "${OUT}/win20.bed"
G1="${OUT}/groups_1.tsv"; awk -F'\t' -v OFS='\t' '{print $1,"G1"}' "${ROOT}/results/staging/real/lyrata_groups.noheader.tsv" > "${G1}"
G12="${IN}/groups_k12.tsv"   # made by scripts/run_a2_real_scaling.sh
# ONLY=<letters> re-times just those stages (e.g. ONLY=f) and replaces their rows in profile.tsv; default: all stages
ONLY="${ONLY:-abcdefg}"
RES="${OUT}/profile.tsv"
if [[ "${ONLY}" == abcdefg || ! -s "${RES}" ]]; then echo -e "stage\tdescription\telapsed_sec" > "${RES}"
else awk -F'\t' -v only="${ONLY}" 'NR==1 || index(only, substr($1,1,1))==0' "${RES}" > "${RES}.tmp" && mv "${RES}.tmp" "${RES}"; fi
t() {  # stage desc cmd...
  s="$1"; d="$2"; shift 2
  [[ "${ONLY}" == *"${s:0:1}"* ]] || return 0
  a=$(perl -MTime::HiRes=time -e 'print time')
  "$@" > /dev/null 2> "${OUT}/${s}.err" || true
  b=$(perl -MTime::HiRes=time -e 'print time')
  echo -e "${s}\t${d}\t$(echo "$b - $a" | bc -l)" | tee -a "${RES}"
}
scan='BEGIN{FS="\t"} /^#/{next} {for(col=10;col<=NF;col++){gtend=index($col,":"); if(!gtend)gtend=length($col)+1; for(c=1;c<gtend;c+=2){al=substr($col,c,1); if(al=="."){m++;break} a[al]++; n++}}} END{print n,m}'
for rep in 1 2 3; do
  t "a${rep}" "tabix only" bash -c "tabix '${VCF}' ${REG} > /dev/null"
  t "b${rep}" "+ field splitting" bash -c "tabix '${VCF}' ${REG} | gawk -F'\t' '{x+=NF} END{print x}'"
  t "c${rep}" "+ genotype scan" bash -c "tabix '${VCF}' ${REG} | gawk '${scan}'"
  t "d${rep}" "piawka pi, 1 group" piawka calc -q -v "${VCF}" -g "${G1}" -b "${OUT}/win20.bed" -s pi
  t "e${rep}" "piawka 5 stats, 1 group" piawka calc -q -v "${VCF}" -g "${G1}" -b "${OUT}/win20.bed" -s pi,theta_w,tajima,lines,miss
  t "f${rep}" "piawka pi, 12 groups" piawka calc -q -v "${VCF}" -g "${G12}" -b "${OUT}/win20.bed" -s pi
  t "g${rep}" "piawka pi, 1 group, UTF-8 locale" env LC_ALL=en_US.UTF-8 piawka calc -q -v "${VCF}" -g "${G1}" -b "${OUT}/win20.bed" -s pi
done
# execution counts (not times) on 2 windows
[[ "${ONLY}" == abcdefg ]] || { echo "profiling done (stages ${ONLY})"; exit 0; }
head -n 2 "${OUT}/win20.bed" > "${OUT}/win2.bed"
( export AWKPATH="${ROOT}/tools/piawka-27e6ae3/include:${AWKPATH:-}"
  gawk --profile="${OUT}/gawk.prof" -E "${ROOT}/tools/piawka-27e6ae3/piawka" calc -q -v "${VCF}" -g "${G1}" -b "${OUT}/win2.bed" -s pi -j 1 > /dev/null 2> "${OUT}/prof.err" ) || true
echo "profiling done"
