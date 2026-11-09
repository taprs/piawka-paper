#!/usr/bin/env bash
# Producer of results/staging/real/northsouth_piawka.bed: the `fst` rows of the default North/South run (run_a3_northsouth.sh) plus
# per-gene `segr2` site counts, obtained by counting the sites of a `piawka --persite -s fst` run (~40 min with 4 jobs).
# Needs results/a3_northsouth/default/piawka.bed (run scripts/run_a3_northsouth.sh first).
# Usage: scripts/run_a3_segr2.sh [jobs=4] [out=results/staging/real/northsouth_piawka.bed]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${ROOT}/tools/piawka-27e6ae3:/Users/ntikhomirov/mambaforge/envs/piawka-paper/bin:${PATH}" TMPDIR=/tmp
J="${1:-4}"; R="${ROOT}/results/staging/real"; OUT="${2:-${R}/northsouth_piawka.bed}"
FST="${ROOT}/results/a3_northsouth/default/piawka.bed"
[[ -s "${FST}" ]] || { echo "missing ${FST}: run scripts/run_a3_northsouth.sh first" >&2; exit 1; }
piawka calc -q -j "${J}" --persite -v "${R}/All_lyrata_allpos.GT.NS.vcf.gz" -g "${R}/northsouth_groups_124.tsv" -b "${R}/genes.bed" -s fst \
  | python3 "${ROOT}/scripts/build_northsouth_piawka_bed.py" "${R}/genes.bed" "${FST}" "${OUT}"
