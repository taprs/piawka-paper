#!/usr/bin/env bash
# Fetches the full allpos VCF (~694 GB!) and its tabix index into $DATA; resumable (curl -C -).
# Skip this and set IN_VCF=<url> to let htslib read the windows remotely by HTTP range requests instead (slower).
set -euo pipefail
source "$(dirname "$0")/config.sh"
for f in "${URL}" "${URL}.tbi"; do
  curl -L -C - --fail --retry 5 --retry-delay 30 -o "${DATA}/$(basename "${f}")" "${f}"
done
ls -l "${DATA}/$(basename "${URL}")"*
