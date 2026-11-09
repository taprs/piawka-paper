#!/usr/bin/env bash
# One window: keep the North/South samples and only GT (plus CHROM/POS/REF/ALT); invariant sites are kept (no -c/-v filters).
# Usage: strip_window.sh <window_number> <contig:start-end>
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "${PARTS}"
out="${PARTS}/$1.vcf.gz"
[[ -s "${out}.ok" ]] && exit 0
bcftools view -S "${SAMPLES_PRESENT}" -r "$2" ${OVL} -Ou "${IN_VCF}" \
  | bcftools annotate -x ID,QUAL,FILTER,INFO,^FORMAT/GT -Oz${ZLEVEL} -o "${out}.tmp"
mv "${out}.tmp" "${out}" && touch "${out}.ok"
