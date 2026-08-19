#!/usr/bin/env bash
# North/South real-data Fst rank comparison (piawka vs pixy) for two gene sets:
# GWAS-SNP-overlapping genes, and GO:2000028 (flowering time) genes.
#
# Starts from pre-staged tool outputs under results/staging/real/ (see
# scripts/analyze_northsouth_fst_ranks.py header for how those were produced
# and what each input is used for). Does not re-run piawka/pixy: the source
# VCF is not staged in this repo (too large) and the pre-staged outputs are
# treated as fixed inputs.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAGE_REAL="${ROOT_DIR}/results/staging/real"

GO_GENES="${STAGE_REAL}/go2000028_athaliana_genes.tsv"
if [[ ! -s "${GO_GENES}" ]]; then
  echo "Fetching GO:2000028 A. thaliana gene list from Ensembl Plants BioMart..."
  bash "${ROOT_DIR}/scripts/fetch_go2000028_athaliana.sh"
fi

ANNOTATION="${STAGE_REAL}/Phytozome/PhytozomeV14/Alyrata/v2.1/annotation/Alyrata_384_v2.1.P14.annotation_info.txt.gz"
if [[ ! -s "${ANNOTATION}" ]]; then
  echo "Missing ${ANNOTATION}" >&2
  echo "Download the Phytozome v14 Alyrata v2.1 annotation_info.txt.gz (JGI Genome Portal login required)" >&2
  echo "and stage it at that path before running this script." >&2
  exit 1
fi

python3 "${ROOT_DIR}/scripts/analyze_northsouth_fst_ranks.py"
Rscript "${ROOT_DIR}/scripts/plot_northsouth_fst_ranks.R"
