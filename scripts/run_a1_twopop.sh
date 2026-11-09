#!/usr/bin/env bash
# A1 (review response): piawka vs pixy on the 16 rescaled two-population SLiM replicates, without and with missing data.
# For each replicate run both tools single-process on one 10 Mbp window and the stats pi, dxy, fst
# (Hudson), Watterson's theta and Tajima's D, for every combination of
#   missingness: none   - the simulated VCF as is
#                masked - every genotype independently set to ./. with a per-site probability drawn from U(0.2, 0.6)
#                         (scripts/add_missing.awk, seed 1000+replicate number, so the masked VCFs are reproducible)
#   variant:     default   - each tool's own default site filtering
#                biallelic - both tools on the same VCF restricted to biallelic SNPs (bcftools view -m2 -M2 -v snps)
# Output: results/a1_twopop16/<Rep>/<none|masked>/<default|biallelic>/{piawka.bed,pixy/}
# Versions: piawka 27e6ae3 (0.9.1 + 2 commits) (tools/piawka-27e6ae3), pixy 2.2.3 (env piawka-paper).
# Usage: scripts/run_a1_twopop.sh [parallel_jobs=6]   (not while timing benchmarks are running)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ENVBIN="${ENVBIN:-/Users/ntikhomirov/mambaforge/envs/piawka-paper/bin}"
export PATH="${ROOT}/tools/piawka-27e6ae3:${ENVBIN}:${PATH}"
export TMPDIR=/tmp
SIM="${ROOT}/results/staging/sim16"
GRP="${SIM}/twopop_groups.tsv"
WIN="${SIM}/onepop_win.bed"
OUT="${ROOT}/results/a1_twopop16"
mkdir -p "${OUT}"

run_one() {
  vcf="$1"; rep="$(basename "${vcf}" .vcf.gz)"; rep="${rep#Demotwopop_scaled10_}"; seed=$((1000 + ${rep#Rep}))
  d="${OUT}/${rep}"; mkdir -p "${d}"
  [[ -s "${d}/done" ]] && return 0
  cp "${GRP}" "${d}/groups.tsv"
  for miss in none masked; do
    src="${vcf}"
    if [[ "${miss}" == masked ]]; then
      src="${d}/masked.vcf.gz"
      bcftools view "${vcf}" | awk -v SEED="${seed}" -f "${ROOT}/scripts/add_missing.awk" | bgzip -c > "${src}"; tabix -f -p vcf "${src}"
    fi
    bcftools view -m2 -M2 -v snps -Oz -o "${d}/biallelic.vcf.gz" "${src}" && tabix -f -p vcf "${d}/biallelic.vcf.gz"
    for variant in default biallelic; do
      in="${src}"; [[ "${variant}" == biallelic ]] && in="${d}/biallelic.vcf.gz"
      o="${d}/${miss}/${variant}"; mkdir -p "${o}"
      piawka calc -q -v "${in}" -g "${d}/groups.tsv" -b "${WIN}" -s pi,dxy,fst,fstwc,theta_w,tajima,lines,miss \
        > "${o}/piawka.bed" 2> "${o}/piawka.err"
      pixy --vcf "${in}" --populations "${d}/groups.tsv" --stats pi dxy fst watterson_theta tajima_d \
        --fst_type hudson --window_size 10000000 --bypass_invariant_check --n_cores 1 \
        --output_folder "${o}/pixy" --output_prefix pixy > "${o}/pixy.log" 2>&1
    done
    rm -f "${d}/biallelic.vcf.gz" "${d}/biallelic.vcf.gz.tbi" "${d}/masked.vcf.gz" "${d}/masked.vcf.gz.tbi"
  done
  date > "${d}/done"
}
export -f run_one
export ROOT SIM GRP WIN OUT
xargs -P "${1:-6}" -I{} bash -c 'run_one {}' < "${SIM}/twopop.txt"
