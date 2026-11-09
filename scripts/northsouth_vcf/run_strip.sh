#!/usr/bin/env bash
# Strips the allpos VCF to North/South samples x GT in parallel windows, then joins and indexes the parts.
# Usage: bash run_strip.sh [parallel_jobs=12]. Needs bcftools, tabix and GNU parallel; safe to re-run after an interruption.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
source ./config.sh
NJOBS="${1:-12}"   # each job loads the index of the huge VCF twice, so keep this modest
[[ -s "${SAMPLES_PRESENT}" ]] || bash check_samples.sh
[[ -s "${DATA}/windows.txt" ]] || bash make_windows.sh
# --resume-failed: on a rerun redo failed windows and run those not started yet
parallel -j "${NJOBS}" --colsep ' ' --joblog "${DATA}/strip.joblog" --resume-failed bash strip_window.sh {1} {2} :::: "${DATA}/windows.txt"
[[ $(ls "${PARTS}"/*.ok | wc -l) -eq $(wc -l < "${DATA}/windows.txt") ]] || { echo "missing parts"; exit 1; }
ls "${PARTS}"/*.vcf.gz | sort > "${DATA}/parts.list"
bcftools concat --naive -f "${DATA}/parts.list" -Oz -o "${OUT_VCF}"
tabix -f -p vcf "${OUT_VCF}"
echo "records in output: $(bcftools index -n "${OUT_VCF}")"; echo "samples: $(bcftools query -l "${OUT_VCF}" | wc -l)"
echo "records in original: $(bcftools index -n "${IN_VCF}")"
ls -l "${OUT_VCF}"
