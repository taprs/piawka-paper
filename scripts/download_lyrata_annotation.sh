#!/usr/bin/env bash
# Download the A. lyrata NT1 reference genome and Helixer gene annotation from https://lyrata.mpipz.mpg.de/downloads
# into results/staging/real/annotation/ (gitignored via results/staging/). Needed by scripts/run_a3_skew_sites.sh
# (bcftools csq). Skips files that are already present (partial downloads use a .part suffix).
# Usage: scripts/download_lyrata_annotation.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BASE="https://lyrata.mpipz.mpg.de/data"
D="${ROOT}/results/staging/real/annotation"; mkdir -p "${D}"
for f in lyrata_NT1.fasta lyrata_NT1.fasta.fai lyrata_genes.gff; do
  [[ -s "${D}/${f}" ]] && { echo "have ${f}"; continue; }
  echo "fetching ${f}"
  curl -sSfL -o "${D}/${f}.part" "${BASE}/${f}" && mv "${D}/${f}.part" "${D}/${f}"
done
ls -l "${D}"
