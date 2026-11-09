#!/usr/bin/env bash
# A3 follow-up: which sites drive the piawka-vs-pixy Fst gap in the genes with non-overlapping jackknife CIs?
# 1. pick the Fig 3 genes (GO:2000028, piawka percentile >= 75, piawka/pixy jackknife CIs not overlapping)
# 2. piawka --persite Hudson Fst (default filtering) on those genes, North vs South (124 individuals)
# 3. bcftools csq consequences for every site in the gene intervals (reference + Helixer GFF from
#    scripts/download_lyrata_annotation.sh)
# Output: results/a3_skew_sites/{genes.bed,persite.bed,sites.csq.vcf}; scripts/analyze_a3_skew_sites.py builds the tables.
# Usage: scripts/run_a3_skew_sites.sh [jobs=4]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${ROOT}/tools/piawka-27e6ae3:/Users/ntikhomirov/mambaforge/envs/piawka-paper/bin:${PATH}" TMPDIR=/tmp
J="${1:-4}"
R="${ROOT}/results/staging/real"; A="${R}/annotation"; O="${ROOT}/results/a3_skew_sites"; mkdir -p "${O}"
VCF="${R}/All_lyrata_allpos.GT.NS.vcf.gz"; G="${R}/northsouth_groups_124.tsv"
for f in lyrata_NT1.fasta lyrata_genes.gff; do
  [[ -s "${A}/${f}" ]] || { echo "missing ${A}/${f}: run scripts/download_lyrata_annotation.sh" >&2; exit 1; }
done
python3 "${ROOT}/scripts/analyze_a3_skew_sites.py" select "${O}/genes.bed"
piawka calc -q -j "${J}" --persite -v "${VCF}" -g "${G}" -b "${O}/genes.bed" -s fst > "${O}/persite.bed" 2> "${O}/persite.err"
# regions are 1-based inclusive: bed start + 1 .. end
awk -F'\t' '{print $1"\t"$2+1"\t"$3}' "${O}/genes.bed" > "${O}/regions.tsv"
# per-population allele counts (INFO AC_<pop>/AN_<pop>), then drop samples for csq (sites-only, every ALT annotated)
bcftools view -R "${O}/regions.tsv" "${VCF}" -Ou \
  | bcftools +fill-tags -Ou -- -S "${G}" -t AN,AC \
  | bcftools view -G -Ou \
  | bcftools csq -f "${A}/lyrata_NT1.fasta" -g "${A}/lyrata_genes.gff" -l -p a --force -Ov -o "${O}/sites.csq.vcf" 2> "${O}/csq.err"
echo "sites: $(grep -vc '^#' "${O}/sites.csq.vcf")"
