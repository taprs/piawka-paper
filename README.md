# piawka paper — benchmarks and figures

Code to reproduce the benchmarks, figures and supplementary tables of:

> Tikhomirov N, Novikova P. *piawka: effective calculator of nucleotide statistics for the thousand-genomes era.*

The manuscript itself is not distributed from this repository.

`piawka` itself lives at <https://github.com/novikovalab/piawka>; this repository only contains the code that
compares it against [`pixy`](https://github.com/ksamuk/pixy). Every analysis here is used in the manuscript.

## What is reproduced

| Manuscript item | Script chain (run in this order) | Aggregated table(s) in `results/tables/` |
| --- | --- | --- |
| **Figure 1** — resource usage on simulated data | `run_resource_sim16.sh` → `analyze_a2_scaling.py` → `make_figures.R` | `a2_scaling.tsv` |
| In-text window-size memory fit | `analyze_a2_window_memory.R` (after `analyze_a2_scaling.py`) | `a2_window_memory_fit.tsv` |
| **Supplementary Figure S1** — resource usage, real and simulated | `run_a2_real_scaling.sh` + `run_resource_sim16.sh` → `analyze_a2_scaling.py` → `plot_a2_scaling_grid.R` | `a2_scaling.tsv` |
| In-text profiling of piawka | `run_a2_profile.sh` → `analyze_a2_scaling.py` | `a2_profile.tsv` |
| **Figure 2** — data gain and π on real data | `run_benchmarks.sh` → `make_figures.R` | `data_gain_pi_vs_site_gain.tsv`, `multiallelic_pi_real_pairs.tsv` |
| In-text CV and correlation tests | `analyze_real_pi_stats_tests.R` | `real_pi_site_gain_stats_tests.tsv` |
| **Supplementary Figure S2** — π ratios, multiallelic vs biallelic | `plot_a4_percent_diff.R` (after `run_benchmarks.sh`) | `multiallelic_pi_ratio_by_population.tsv` |
| In-text multiallelic π uplift on simulated data | `run_resource_sim16.sh` → `run_benchmarks.sh` | `multiallelic_pi_simulated_pairs.tsv` |
| In-text missing-data consistency, **Supplementary Table 2** | `run_a1_twopop.sh`, `analyze_a1_consistency.py`, `collect_fst_tables.py` | `a1_twopop_{per_rep,summary}.tsv` |
| **Figure 3**, **Supplementary Figure S3**, **Supplementary Tables 3–4** | Gene-wise Fst chain below | `northsouth_fst_ranks_*`, `northsouth_jackknife_ci_*` |
| In-text identical-sites check | `run_a3_northsouth.sh` → `analyze_a3_identical_sites.py` | `northsouth_identical_sites.tsv` |
| **Supplementary Table 5** — sites driving the Fst gap | `download_lyrata_annotation.sh`, `run_a3_skew_sites.sh`, `analyze_a3_skew_sites.py`, `collect_fst_tables.py` | `northsouth_skew_*.tsv` |

All aggregated tables are committed, so **figures can be rebuilt without re-running any benchmark**
(see [Rebuilding the figures](#rebuilding-the-figures-from-the-committed-tables)). Re-running the benchmarks needs the input
data ([Input data](#input-data)). Raw run outputs under `results/` (other than `results/tables/`) are git-ignored.

## Setup

```bash
conda env create -f env/environment.yml        # env name: piawka-paper
# cvequality is CRAN-only; needed by scripts/analyze_real_pi_stats_tests.R
conda run -n piawka-paper Rscript -e 'install.packages("cvequality", repos="https://cloud.r-project.org")'
```

- **piawka** is pinned to commit `27e6ae3` (tag 0.9.1 plus two commits that make Watterson's θ and Tajima's D aggregate sites with
  unequal numbers of called genotypes the way pixy does; not on bioconda, whose latest is 0.9.0). The scripts take it from
  `tools/piawka-27e6ae3` (git-ignored):

  ```bash
  mkdir -p tools && git clone https://github.com/novikovalab/piawka.git tools/piawka-27e6ae3
  git -C tools/piawka-27e6ae3 checkout 27e6ae3
  ```

- **pixy 2.2.3** is an untagged `master` build, pinned by commit (`21c7c33`) in `env/environment.yml`.
- Scripts that put the tools on `PATH` themselves read the environment's `bin/` from `ENVBIN`
  (default `/Users/ntikhomirov/mambaforge/envs/piawka-paper/bin`); set it to your own. `run_benchmarks.sh` does not:
  run it with `PATH="$PWD/tools/piawka-27e6ae3:<env>/bin:$PATH"`.
- Benchmarks ran on a MacBook Pro M1 with 32 GB RAM. The timing wrappers use BSD `/usr/bin/time -l` and `ps`, so they are
  written for macOS; aggregation and plotting are portable. Timing runs are sequential and need an otherwise idle machine.

## Rebuilding the figures from the committed tables

```bash
conda run -n piawka-paper Rscript scripts/make_figures.R                      # figures/main/f1.png, f2.png
conda run -n piawka-paper Rscript scripts/plot_a2_scaling_grid.R              # figures/supp/S1_resource_usage.png
conda run -n piawka-paper Rscript scripts/plot_a4_percent_diff.R              # figures/supp/S2_pi_ratios.png
conda run -n piawka-paper Rscript scripts/plot_northsouth_fst_ranks_jackknife.R   # figures/main/f3.png, figures/supp/S3_flowering_time_fst.png
conda run -n piawka-paper Rscript scripts/analyze_a2_window_memory.R          # window-memory fit
conda run -n piawka-paper Rscript scripts/analyze_real_pi_stats_tests.R       # in-text tests
conda run -n piawka-paper python scripts/collect_fst_tables.py                # Supplementary Tables 1-5 workbook
```

## Reproducing the analyses

### Resource usage (Figure 1, Supplementary Figure S1)

```bash
bash scripts/run_resource_sim16.sh       # 16 neutral replicates: window sweep, 1/2/4/8 processes, 10-100 samples
bash scripts/run_a2_real_scaling.sh      # real data: 25-910 samples at 1 and 8 processes; 1-12 populations
bash scripts/run_a2_profile.sh           # piawka stage ablation on 200 kb of the real data
python scripts/analyze_a2_scaling.py     # -> a2_scaling.tsv, a2_profile.tsv  (needs pandas)
```

`run_resource_sim16.sh` writes `results/resource_sim16/<Rep>/<tag>/`: `win<bp>_p1` (10 Mbp analysed as windows of
0.625–10 Mbp, one process), `win1250000_p{1,2,4,8}` (eight 1.25 Mbp windows dispatched to 1/2/4/8 processes, the process
series of Figure 1B; neither tool splits a single window between processes) and `n<k>_p{1,8}` (nested random subsets of
10–75 samples with sites that became monomorphic dropped; 100 samples are the `win1250000` runs).

### Data gain, multiallelic π (Figure 2, S2, simulated in-text claim)

```bash
PATH="$PWD/tools/piawka-27e6ae3:<env>/bin:$PATH" bash scripts/run_benchmarks.sh
```

Runs, then aggregates with `analyze_benchmarks.py`: π in 100 windows of the real *A. lyrata* 1 Mbp slice with both tools
(data gain), the same with `piawka --mult` / `pixy --include_multiallelic_snps`, and `piawka --mult` on the 16 neutral
replicates compared with the biallelic runs of `run_resource_sim16.sh` (run that first). Skips runs that already exist.

### Missing-data consistency (Supplementary Table 2)

```bash
bash scripts/run_a1_twopop.sh                      # 16 two-population replicates x {none, masked} x {default, biallelic}
python scripts/analyze_a1_consistency.py           # -> a1_twopop_{per_rep,summary}.tsv
```

Masking is `add_missing.awk` (every genotype set to `./.` with a per-site probability drawn from U(0.2, 0.6), seed 1000 +
replicate number).

### Gene-wise North/South Fst (Figure 3, S3, Supplementary Tables 3–5)

Compares piawka and pixy per-gene Hudson's Fst **percentile ranks** between northern (117) and southern (7) individuals of
the East Siberian *A. lyrata* lineage, for genes overlapping environmental-GWAS SNPs and for genes whose *A. thaliana*
homolog carries GO:2000028 ("regulation of photoperiodism, flowering").

0. `python scripts/make_northsouth_groups.py` — the 124 East Siberian diploid individuals from the public metadata table (`Admixture population` = East Siberia, ploidy 2), split into "South" (latitude < 57, 7 individuals) and "North" (117) → `results/staging/real/northsouth_groups_124.tsv` and `scripts/northsouth_vcf/samples_ns.txt`.
1. Build the GT-only North/South VCF from the full-genome VCF (`scripts/northsouth_vcf/`, see its README).
2. `bash scripts/run_a3_northsouth.sh` — piawka and pixy on it, default filtering and both tools on the same biallelic
   SNPs → `results/a3_northsouth/`. `python scripts/analyze_a3_identical_sites.py` → `northsouth_identical_sites.tsv`.
   `bash scripts/run_a3_segr2.sh` — `northsouth_piawka.bed` (the `fst` rows of the default run plus per-gene `segr2` site counts from a `--persite` run).
3. `bash scripts/run_northsouth_fst_comparison.sh` — rank tables `northsouth_fst_ranks_{gwas_snp,go2000028}_genes.tsv`,
   from the staged `northsouth_piawka.bed` / `northsouth_pixy_fst.txt` (see input table). It bridges *A. lyrata* gene IDs
   to *A. thaliana* orthologs with the login-gated JGI Phytozome `Alyrata_384_v2.1.P14.annotation_info.txt.gz` (stage it
   at `results/staging/real/Phytozome/PhytozomeV14/Alyrata/v2.1/annotation/`); the GO:2000028 list is fetched by
   `fetch_go2000028_athaliana.sh`. The docstring of `analyze_northsouth_fst_ranks.py` explains the rank definitions.
4. `python scripts/analyze_northsouth_diagnostics.py` — genome-wide diagnostics and `northsouth_fst_all_genes.tsv`.
5. `bash scripts/run_a3_jackknife.sh` — leave-one-South-out folds (7 folds, 117 North vs 6 South); then
   `python scripts/analyze_a3_jackknife.py` → `northsouth_jackknife_ci_*_genes.tsv`.
6. `Rscript scripts/plot_northsouth_fst_ranks_jackknife.R` → `figures/main/f3.png`, `figures/supp/S3_flowering_time_fst.png`.
7. Supplementary Table 5: `bash scripts/download_lyrata_annotation.sh` (reference + Helixer GFF),
   `bash scripts/run_a3_skew_sites.sh`, `python scripts/analyze_a3_skew_sites.py`.
8. `python scripts/collect_fst_tables.py` → `results/tables/piawka_supplementary_data.xlsx` (Supplementary Tables 1–5 and a
   column legend; needs openpyxl).

## Input data

Inputs and intermediates are not committed (multiple GB); `results/staging/` and the raw run directories are git-ignored.

| Path under `results/staging/` | What it is |
| --- | --- |
| `sim16/onepop/<Rep>.vcf.gz`, `sim16/onepop.txt` | the 16 one-population VCFs (`Demoonepop_scaled10_Rep0–15`, 100 samples, ~460k records) and the list of their paths (relative to the project root) |
| `sim16/onepop_groups.tsv`, `sim16/onepop_win.bed` | sample → population map (all `POP1`) and one 10 Mbp window (`1 0 10000000`) |
| `sim16/twopop/Demotwopop_scaled10_Rep{0–15}.vcf.gz`, `sim16/twopop.txt`, `sim16/twopop_groups.tsv` | the 16 two-population VCFs (100 samples), their list, and the sample map (`i0–i49` = p1, `i50–i99` = p2) |
| `real/All_lyrata_final_allpos_1mbp.vcf.gz` (+ `.tbi`) | real *A. lyrata* VCF, `scaffold_1:1000001-2000000` |
| `real/lyrata_win.bed`, `real/lyrata_groups.tsv` | 100 windows of 10 kbp; ADMIXTURE population labels for the 910 *A. lyrata* samples (12 labels, identical to the public metadata table). Figure 2 drops the `MIX*`/`admixed` labels (77 samples), leaving 833 samples in eight lineages |
| `real/All_lyrata_allpos.GT.NS.vcf.gz` | GT-only VCF of the 124 North/South individuals (step 1 above); `real/northsouth_groups_124.tsv` is its group table |
| `real/genes.bed`, `real/supplementary_data_4_Lyrata_GWAS_SNPs.tsv`, `real/northsouth_piawka.bed`, `real/northsouth_pixy_fst.txt` | inputs of the rank tables. The pixy file equals the `default` run of `run_a3_northsouth.sh`; the piawka file has the same `fst` rows plus per-gene `segr2` rows (segregating sites usable by the pair; they feed the site counts in Supplementary Tables 3–4), rebuilt by `run_a3_segr2.sh` from `piawka --persite` counts. The submitted file also had `lines2` rows from an older build; nothing reads them |
| `real/annotation/` | reference genome and Helixer GFF (`download_lyrata_annotation.sh`) |

The real dataset is from Glushkevich et al. (2026, bioRxiv 2026.03.02.709016), accessed via
<https://www.arabidopsislyrata.org/>.

## Simulated data (`slim/`)

The VCFs are not committed, but their SLiM (3.7.1) recipes are: one directory per replicate under `slim/`, each with the
`infile_nuc.txt` that produced that replicate's VCF (`slim infile_nuc.txt` writes `out.vcf`; bgzip and tabix it).

- `Demoonepop_scaled10_Rep0–15` — one constant-size population, 100 diploids sampled at generation 5,500.
- `Demotwopop_scaled10_Rep0–15` — two populations that split completely after generation 5,000, 50 diploids sampled from each.

Both follow Zeitler & Gilbert (2024, *GBE* 16:evae139) with a nucleotide-based model (Jukes–Cantor substitutions,
`simplifyNucleotides=T` on output) and are rescaled tenfold to run faster: 1,000 individuals, 5,500 generations, mutation
rate 7 × 10⁻⁷ and recombination rate 1 × 10⁻⁷ per site per generation, one 10 Mbp chromosome. Every recipe carries an
explicit `setSeed(...)`.

## Repository layout

```
env/environment.yml     pinned conda environment
slim/                   SLiM recipes for the 16 + 16 simulated replicates
scripts/                run_* (benchmarks), analyze_* (aggregation), plot_* / make_figures.R (figures)
results/tables/         aggregated TSVs and the supplementary workbook (committed)
figures/main/           f1.png, f2.png, f3.png
figures/supp/           S1_resource_usage.png, S2_pi_ratios.png, S3_flowering_time_fst.png
```

| Script | Role |
| --- | --- |
| `run_piawka_benchmark.sh`, `run_pixy_benchmark.sh` | single-run wrappers recording wall time, CPU time and peak RSS over the whole process group |
| `run_resource_sim16.sh`, `run_a2_real_scaling.sh`, `run_a2_profile.sh` | resource usage on simulated data, on real data, and the profiling ablation |
| `analyze_a2_scaling.py`, `analyze_a2_window_memory.R`, `plot_a2_scaling_grid.R` | aggregate the resource runs; memory-vs-window fit; Supplementary Figure S1 |
| `run_benchmarks.sh`, `analyze_benchmarks.py` | data gain and multiallelic π (real and simulated), and their aggregation |
| `make_figures.R` | Figures 1 and 2 |
| `analyze_real_pi_stats_tests.R`, `plot_a4_percent_diff.R` | in-text tests on Figure 2 data; Supplementary Figure S2 |
| `run_a1_twopop.sh`, `add_missing.awk`, `analyze_a1_consistency.py` | missing-data consistency on the two-population replicates |
| `northsouth_vcf/`, `run_a3_northsouth.sh`, `analyze_a3_identical_sites.py` | North/South GT-only VCF; default vs identical-site Fst runs |
| `make_northsouth_groups.py`, `run_a3_segr2.sh`, `build_northsouth_piawka_bed.py` | North/South sample groups; `northsouth_piawka.bed` |
| `run_northsouth_fst_comparison.sh`, `analyze_northsouth_fst_ranks.py`, `fetch_go2000028_athaliana.sh` | gene-wise Fst rank tables |
| `analyze_northsouth_diagnostics.py`, `run_a3_jackknife.sh`, `analyze_a3_jackknife.py` | genome-wide diagnostics; leave-one-South-out confidence intervals |
| `plot_northsouth_fst_ranks_jackknife.R` | Figure 3 and Supplementary Figure S3 |
| `download_lyrata_annotation.sh`, `run_a3_skew_sites.sh`, `analyze_a3_skew_sites.py` | sites driving the Fst gap (Supplementary Table 5) |
| `collect_fst_tables.py` | Supplementary Tables 1–5 workbook |
