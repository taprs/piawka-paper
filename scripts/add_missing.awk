#!/usr/bin/env -S mawk -f

# This script replaces a share of samples' genotypes (between MIN and MAX) at each site by missing data. For testing purposes
# Can be sped up with GNU parallel: zcat $vcf | parallel --pipe --keep-order --block 10M add_missing.awk | bgzip > output.vcf.gz

# Arguments: MIN and MAX between 0 and 1 (default 0.2 and 0.6); SEED (integer) for reproducible masking

BEGIN{ srand(SEED ? SEED : srand()); OFS="\t"
  if ( !MIN ) { MIN=0.2 }
  if ( !MAX ) { MAX=0.6 } 
}
/^[^#]/ {
		for (i=10; i<=NF; i++) {
      # Probability threshold between MAX and MIN
      prob=MIN+rand()*(MAX-MIN)

      if ( rand() < prob ) {
        gt = "./."
        #gt = "."
        #for (j=1; j < match($i, ":") / 2; j++) { gt = gt"/." }
        $i = gt #substr($i, match($i, ":"))
        }
			}
}1

