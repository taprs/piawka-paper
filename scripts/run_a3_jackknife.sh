#!/usr/bin/env bash
# Leave-one-South-out jackknife (7 South individuals -> 7 folds, 117 North vs 6 South; results/a3_jackknife/groups/jkNN.tsv):
# gene-wise FST with piawka 27e6ae3 (0.9.1 + 2 commits) and pixy 2.2.3 under default site filtering, for the Fig 3 CI variant.
# Usage: scripts/run_a3_jackknife.sh [parallel_folds=2] [jobs_per_run=2]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${ROOT}/tools/piawka-27e6ae3:/Users/ntikhomirov/mambaforge/envs/piawka-paper/bin:${PATH}" TMPDIR=/tmp
export ROOT
R="${ROOT}/results/staging/real"; export VCF="${R}/All_lyrata_allpos.GT.NS.vcf.gz" B="${R}/genes.bed" J="${2:-2}"
run_one() {
  g="$1"; n="$(basename "${g}" .tsv)"; d="${ROOT}/results/a3_jackknife/${n}"; mkdir -p "${d}/pixy"
  [[ -s "${d}/done" ]] && return 0
  piawka calc -q -j "${J}" -v "${VCF}" -g "${g}" -b "${B}" -s fst > "${d}/piawka.bed" 2> "${d}/piawka.err"
  pixy --vcf "${VCF}" --populations "${g}" --stats fst --fst_type hudson --bed_file "${B}" --bypass_invariant_check \
    --n_cores "${J}" --output_folder "${d}/pixy" --output_prefix pixy > "${d}/pixy.log" 2>&1
  date > "${d}/done"
}
export -f run_one
# fold jkNN = the 124-individual group table without its NN-th South individual (file order)
GD="${ROOT}/results/a3_jackknife/groups"; mkdir -p "${GD}"
n_south="$(awk -F'\t' '$2=="South"' "${R}/northsouth_groups_124.tsv" | wc -l | tr -d ' ')"
for i in $(seq 1 "${n_south}"); do
  f="${GD}/jk$(printf '%02d' "${i}").tsv"
  [[ -s "${f}" ]] || awk -F'\t' -v k="${i}" '$2=="South"{n++; if (n==k) next} {print}' "${R}/northsouth_groups_124.tsv" > "${f}"
done
ls "${ROOT}"/results/a3_jackknife/groups/jk*.tsv | xargs -P "${1:-2}" -I{} bash -c 'set -euo pipefail; run_one {}'
