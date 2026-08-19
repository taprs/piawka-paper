#!/usr/bin/env bash
# Fetches the A. thaliana gene list annotated with GO:2000028 (regulation of
# photoperiodism, flowering) from the Ensembl Plants BioMart REST API, and
# caches it as a staged input. Network access required; re-run only if the
# GO annotation needs refreshing.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT_DIR}/results/staging/real/go2000028_athaliana_genes.tsv"
GO_ID="${1:-GO:2000028}"

QUERY_XML=$(cat <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="plants_mart" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
  <Dataset name="athaliana_eg_gene" interface="default">
    <Filter name="go" value="${GO_ID}"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
  </Dataset>
</Query>
EOF
)

mkdir -p "$(dirname "${OUT}")"
curl -sf --data-urlencode "query=${QUERY_XML}" "https://plants.ensembl.org/biomart/martservice" > "${OUT}"

n=$(tail -n +2 "${OUT}" | wc -l | tr -d ' ')
echo "Wrote ${n} ${GO_ID} A. thaliana genes to ${OUT}"
