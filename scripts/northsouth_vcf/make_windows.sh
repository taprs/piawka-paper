#!/usr/bin/env bash
# Writes windows.txt: "<5-digit number> <contig:start-end>" per window, in the contig order of the VCF index.
set -euo pipefail
source "$(dirname "$0")/config.sh"
bcftools index -s "${IN_VCF}" | awk -v W="${WINDOW}" '$2 ~ /^[0-9]+$/ {for(s=1;s<=$2;s+=W){e=s+W-1; if(e>$2)e=$2; print $1":"s"-"e}}' \
  | nl -n rz -w5 -s' ' > "${DATA}/windows.txt"
echo "$(wc -l < "${DATA}/windows.txt") windows written to ${DATA}/windows.txt"
