# North/South GT-only VCF (analysis A3)

Reduces the 694 GB `All_lyrata_final_allpos.vcf.gz` (1018 samples, all sites) to the 124 North/South individuals
(`samples_ns.txt`) with GT only, invariant sites kept, for the identical-site piawka/pixy FST comparison.

    bash scripts/northsouth_vcf/download.sh        # optional, ~694 GB; or set IN_VCF=<url> to read remotely
    bash scripts/northsouth_vcf/test_one_window.sh # must print OK
    bash scripts/northsouth_vcf/run_strip.sh 12    # parallel windows, concat, tabix

Output: `results/staging/real/All_lyrata_allpos.GT.NS.vcf.gz` (+`.tbi`), ~3.2 GB (parts are bgzip level 9;
level 1 gave 4.9 GB). Override `DATA`, `IN_VCF`, `OUT_VCF`, `URL` via environment (see `config.sh`).
