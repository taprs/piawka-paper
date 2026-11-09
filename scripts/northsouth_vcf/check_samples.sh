#!/usr/bin/env bash
# Compares samples_ns.txt with the VCF header; writes samples_present.txt and samples_missing.txt and reports.
set -euo pipefail
source "$(dirname "$0")/config.sh"
bcftools query -l "${IN_VCF}" | sort > "${DATA}/vcf_samples.txt"
sort "${SAMPLES}" > "${DATA}/.want.txt"
comm -12 "${DATA}/.want.txt" "${DATA}/vcf_samples.txt" > "${SAMPLES_PRESENT}"
comm -23 "${DATA}/.want.txt" "${DATA}/vcf_samples.txt" > "${DATA}/samples_missing.txt"
echo "wanted $(wc -l < "${SAMPLES}"), in VCF $(wc -l < "${SAMPLES_PRESENT}"), missing $(wc -l < "${DATA}/samples_missing.txt"):"
cat "${DATA}/samples_missing.txt"
