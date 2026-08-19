# piawka paper — benchmarks and figures

Code to reproduce the benchmarks and figures in:

> Tikhomirov N, Novikova P. *piawka: effective calculator of nucleotide statistics for the thousand-genomes era.*

The manuscript as submitted is in [`manuscript/piawka_paper.pdf`](manuscript/piawka_paper.pdf).

`piawka` itself lives at <https://github.com/novikovalab/piawka> — this repository only contains the
benchmarking code that compares it against [`pixy`](https://github.com/ksamuk/pixy).

## What is reproduced

| Paper item | Produced by | From table(s) |
| --- | --- | --- |
| **Figure 1** — resource usage on simulated data | `scripts/make_figures.R` | `resource_usage_1thread_synwin_vs_threads8.tsv` (1A), `resource_usage_timing.tsv` (1B) |
| **Figure 2** — data gain and π on real data | `scripts/make_figures.R` | `data_gain_pi_vs_site_gain.tsv`, `multiallelic_pi_real_pairs.tsv` |
| **Figure 3** — North/South gene-wise F<sub>ST</sub> | `scripts/plot_northsouth_fst_ranks.R` | `northsouth_fst_ranks_{gwas_snp,go2000028}_genes.tsv` |
| **Supplementary Tables 1 & 2** | `scripts/analyze_northsouth_fst_ranks.py` | the same two `northsouth_fst_ranks_*.tsv` files |
| In-text statistics (Spearman ρ, CV, Feltz & Miller test) | `scripts/analyze_real_pi_stats_tests.R` | `real_pi_site_gain_stats_tests.tsv` |
| In-text claim: multiallelic π uplift on simulated data | `scripts/analyze_benchmarks.py` | `multiallelic_pi_simulated_pairs.tsv` |

All aggregated tables are committed under `results/tables/`, so **the figures can be rebuilt without
re-running any benchmark**. Re-running the benchmarks themselves additionally requires the input data
(see [Input data](#input-data)).

## Setup

```bash
conda env create -f env/environment.yml
# cvequality is CRAN-only; needed by scripts/analyze_real_pi_stats_tests.R
conda run -n piawka-paper Rscript -e 'install.packages("cvequality", repos="https://cloud.r-project.org")'
```

The environment pins the exact tool versions benchmarked in the paper: **piawka 0.9.0** (bioconda) and
**pixy 2.0.0.dev2**. The latter is an untagged `master` build, so it is pinned by commit
(`fe3c9058`); the nearest tags around it are `2.0.0.beta14`–`2.0.0.beta17`.

Benchmarks in the paper were run on a MacBook Pro M1 with 32 GB RAM. The timing wrappers use BSD
`/usr/bin/time -l` and `ps`, so they are written for macOS; the analysis and plotting steps are portable.

## Reproducing the figures

With the committed tables, this is all that is needed:

```bash
conda run -n piawka-paper Rscript scripts/make_figures.R              # -> figures/main/f1.png, f2.png
conda run -n piawka-paper Rscript scripts/plot_northsouth_fst_ranks.R # -> figures/main/f3.png
conda run -n piawka-paper Rscript scripts/analyze_real_pi_stats_tests.R
```

## Reproducing the benchmarks

### 1. Simulated + real-data benchmarks (Figures 1 and 2)

```bash
conda run -n piawka-paper bash scripts/run_benchmarks.sh
```

This runs, then aggregates into `results/tables/` via `scripts/analyze_benchmarks.py`:

1. **Resource usage** (Figure 1B) — π on each simulated VCF with piawka (`-s pi,lines,miss`) and pixy
   (`--stats pi`), at 1, 2, 4 and 8 parallel processes, over one 10 Mbp window.
2. **Window-size comparison** (Figure 1A) — the same, at a single process, over one 10 Mbp window
   (`syn_win.bed`) versus eight equal 1.25 Mbp windows (`syn_win_threads8.bed`). Run on its own, not
   concurrently with other stages, so the memory measurements are not contended.
3. **Data gain** (Figure 2A, 2B) — π in 10 kbp windows on the real *A. lyrata* VCF with both tools.
4. **Multiallelic** (Figure 2C) — the same real windows plus all simulated VCFs, with `piawka --mult`
   and `pixy --include_multiallelic_snps`.

One caveat on exact reproduction: the committed `resource_usage_timing.tsv` comes from runs that also
computed `theta_w` (`-s pi,lines,miss,theta_w`), whereas the script now uses the command reported in the
paper (`-s pi,lines,miss`). A re-run will therefore give slightly faster piawka timings than the committed
table. This does not affect any conclusion.

To re-aggregate only, once raw outputs exist under `results/`:

```bash
conda run -n piawka-paper python scripts/analyze_benchmarks.py
```

### 2. North/South F<sub>ST</sub> gene-rank comparison (Figure 3, Supplementary Tables 1–2)

An independent analysis, not part of `run_benchmarks.sh`. It compares piawka vs pixy per-gene Hudson's
F<sub>ST</sub> **percentile ranks** between northern and southern individuals of the East Siberian
*A. lyrata* lineage, for two gene sets: genes overlapping environmental-GWAS SNPs, and genes whose
*A. thaliana* homolog carries GO:2000028 ("regulation of photoperiodism, flowering").

```bash
conda run -n piawka-paper bash scripts/run_northsouth_fst_comparison.sh
```

It starts from pre-staged tool outputs (`results/staging/real/northsouth_piawka.bed`,
`northsouth_pixy_fst.txt`) and does **not** re-run piawka/pixy — the source VCF is too large to stage here.

Bridging *A. lyrata* gene IDs (`AL#G#####`) to *A. thaliana* orthologs requires the JGI Phytozome
annotation file `Alyrata_384_v2.1.P14.annotation_info.txt.gz` (login-gated; stage it manually at
`results/staging/real/Phytozome/PhytozomeV14/Alyrata/v2.1/annotation/`). The GO:2000028 gene list is
fetched automatically from Ensembl Plants BioMart by `scripts/fetch_go2000028_athaliana.sh`. See the
docstring of `scripts/analyze_northsouth_fst_ranks.py` for input provenance, rank/percentile definitions,
and why positional overlap against Ensembl's gene models was rejected in favour of the Phytozome crosswalk.

## Input data

Input and intermediate data are **not** committed (they are multiple GB); `results/staging/` and the
raw per-run output directories are git-ignored. To re-run the benchmarks you need to stage:

| Path under `results/staging/` | What it is |
| --- | --- |
| `simulated/onepop.txt`, `simulated/bneck.txt` | lists of paths to the simulated VCFs (64 replicates each of the neutral and long-bottleneck scenarios) |
| `simulated/syn_win.bed`, `simulated/syn_win_threads8.bed` | window definitions (one 10 Mbp window; eight 1.25 Mbp windows) |
| `simulated/syn_groups.tsv` | sample→population map for the simulated data (all one population) |
| `real/All_lyrata_final_allpos_1mbp.vcf.gz` (+ `.tbi`) | real *A. lyrata* VCF, `scaffold_1:1000001-2000000` |
| `real/lyrata_win.bed`, `real/lyrata_groups.tsv` | 100 equal windows; ADMIXTURE population labels for the 815 analyzed samples |
| `real/northsouth_piawka.bed`, `real/northsouth_pixy_fst.txt`, `real/genes.bed`, `real/supplementary_data_4_Lyrata_GWAS_SNPs.tsv` | inputs to the Figure 3 analysis |

Simulations follow Zeitler & Gilbert (2024, *GBE* 16:evae139) with a nucleotide model, a Jukes–Cantor
substitution matrix, and `simplifyNucleotides` enabled. The real dataset is from Glushkevich et al.
(2026, bioRxiv 2026.03.02.709016), accessed via <https://www.arabidopsislyrata.org/>.

## Repository layout

```
env/environment.yml   pinned conda environment
scripts/              benchmark wrappers, aggregation, plotting
results/tables/       aggregated, analysis-ready TSVs (committed)
figures/main/         f1.png, f2.png, f3.png
manuscript/           submitted PDF
```

| Script | Role |
| --- | --- |
| `run_benchmarks.sh` | orchestrates the simulated + real benchmarks, then aggregates |
| `run_piawka_benchmark.sh`, `run_pixy_benchmark.sh` | single-run wrappers recording wall time, CPU time and peak RSS across the whole process group |
| `analyze_benchmarks.py` | aggregates raw run outputs into `results/tables/` (Python stdlib only) |
| `make_figures.R` | Figures 1 and 2 |
| `run_northsouth_fst_comparison.sh` | driver for the Figure 3 analysis |
| `fetch_go2000028_athaliana.sh` | fetches the GO:2000028 gene list from Ensembl Plants BioMart |
| `analyze_northsouth_fst_ranks.py` | builds the gene-wise F<sub>ST</sub> rank tables (Supplementary Tables 1–2) |
| `plot_northsouth_fst_ranks.R` | Figure 3 |
| `analyze_real_pi_stats_tests.R` | in-text correlation and CV-equality tests |
