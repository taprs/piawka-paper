#!/usr/bin/env bash
# Run this first. Strips window 1 and checks it against the original (same records, same genotypes).
set -euo pipefail
source "$(dirname "$0")/config.sh"
bash "$(dirname "$0")/check_samples.sh"
[[ -s "${DATA}/windows.txt" ]] || bash "$(dirname "$0")/make_windows.sh"
read -r n reg < <(sed -n 1p "${DATA}/windows.txt")
echo "bcftools: $(bcftools --version | head -1); overlap flag: '${OVL}'; window ${n} ${reg}"
bash "$(dirname "$0")/strip_window.sh" "${n}" "${reg}"
echo "samples in part: $(bcftools query -l "${PARTS}/${n}.vcf.gz" | wc -l) (expect $(wc -l < "${SAMPLES_PRESENT}"))"
a=$(bcftools query -S "${SAMPLES_PRESENT}" -r "${reg}" ${OVL} -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "${IN_VCF}" | cksum)
b=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "${PARTS}/${n}.vcf.gz" | cksum)
[[ "${a}" == "${b}" ]] && echo "OK: records and genotypes identical to the original" || { echo "MISMATCH"; exit 1; }
ls -l "${PARTS}/${n}.vcf.gz"
