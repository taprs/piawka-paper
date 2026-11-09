#!/usr/bin/env bash
# A3 (review response): gene-wise North/South FST, piawka 27e6ae3 (0.9.1 + 2 commits) vs pixy 2.2.3, on the stripped GT-only VCF
# (scripts/northsouth_vcf; 124 individuals incl. BAM_12.1-1), two variants:
#   default    - each tool's own site filtering on the all-sites VCF
#   biallelic  - both tools on the same VCF restricted to biallelic SNPs (bcftools view -m2 -M2 -v snps)
# Output: results/a3_northsouth/<variant>/{piawka.bed,pixy/pixy_fst.txt,*.time}
# Usage: scripts/run_a3_northsouth.sh [jobs=4]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${ROOT}/tools/piawka-27e6ae3:/Users/ntikhomirov/mambaforge/envs/piawka-paper/bin:${PATH}" TMPDIR=/tmp
J="${1:-4}"
R="${ROOT}/results/staging/real"; O="${ROOT}/results/a3_northsouth"; mkdir -p "${O}"
VCF="${R}/All_lyrata_allpos.GT.NS.vcf.gz"; G="${R}/northsouth_groups_124.tsv"; B="${R}/genes.bed"
if [[ ! -s "${O}/biallelic.vcf.gz.tbi" ]]; then
  bcftools view --threads "${J}" -m2 -M2 -v snps -Oz -o "${O}/biallelic.vcf.gz" "${VCF}" && tabix -f -p vcf "${O}/biallelic.vcf.gz"
fi
for variant in default biallelic; do
  in="${VCF}"; [[ "${variant}" == biallelic ]] && in="${O}/biallelic.vcf.gz"
  d="${O}/${variant}"; mkdir -p "${d}/pixy"
  [[ -s "${d}/piawka.done" ]] || { /usr/bin/time -l -o "${d}/piawka.time" piawka calc -q -j "${J}" -v "${in}" -g "${G}" -b "${B}" -s fst,lines > "${d}/piawka.bed" 2> "${d}/piawka.err" && touch "${d}/piawka.done"; }
  [[ -s "${d}/pixy.done" ]] || { /usr/bin/time -l -o "${d}/pixy.time" pixy --vcf "${in}" --populations "${G}" --stats fst --fst_type hudson \
      --bed_file "${B}" --bypass_invariant_check --n_cores "${J}" --output_folder "${d}/pixy" --output_prefix pixy > "${d}/pixy.log" 2>&1 && touch "${d}/pixy.done"; }
done
