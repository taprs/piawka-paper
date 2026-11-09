# Shared settings for the North/South GT-only extraction (review analysis A3). Edit here or override via environment.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
URL="${URL:-https://lyrata.mpipz.mpg.de/downloads/All_lyrata_final_allpos.vcf.gz}"
DATA="${DATA:-${ROOT}/results/staging/real/ns_strip}"      # scratch/work directory (raw download ~694 GB, parts, logs)
IN_VCF="${IN_VCF:-${DATA}/$(basename "${URL}")}"
SAMPLES="$(dirname "${BASH_SOURCE[0]}")/samples_ns.txt"    # 117 North + 7 South (the VCF calls BAM12.1-1 "BAM_12.1-1")
WINDOW=5000000                                             # window size (bp) for the parallel split
OUT_VCF="${OUT_VCF:-${ROOT}/results/staging/real/All_lyrata_allpos.GT.NS.vcf.gz}"
ZLEVEL=9                                                   # bgzip level of the parts (kept by the naive concat)
PARTS="${DATA}/parts"
SAMPLES_PRESENT="${DATA}/samples_present.txt"              # written by check_samples.sh
mkdir -p "${DATA}"
# bcftools >= 1.15 can restrict a record to the window containing its start position; older versions cannot
OVL=""
if bcftools --version | head -1 | awk '{split($2,v,"."); exit !(v[1]>1 || (v[1]==1 && v[2]>=15))}'; then OVL="--regions-overlap pos"; fi
# conda bcftools links OpenBLAS, which starts one thread per core in EVERY process; with many parallel jobs that exhausts
# the per-user thread limit ("pthread_create failed"). Keep each process single-threaded.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GOTO_NUM_THREADS=1 MKL_NUM_THREADS=1
