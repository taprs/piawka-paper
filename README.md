# piawka_paper benchmark workspace

This repository contains the benchmark and figure-generation workflow for evaluating `piawka` against `pixy` in a methods-paper setting.

## What this benchmark is for

The benchmark is designed to answer four practical questions:

1. **Performance:** How do runtime and memory use compare between `piawka` and `pixy` on the same simulated inputs?
2. **Missing-data stability:** How stable are `theta_watterson` and selection-proxy differences (`pi - theta_w`, `pi - theta_low`) as missingness increases?
3. **Data gain on real data:** How many additional callable sites/windows does `piawka` retain relative to `pixy`, and how does that relate to estimated `pi`?
4. **Multiallelic behavior:** How do `piawka` biallelic and `piawka --mult` estimates compare in simulated and real datasets?

## What is performed

`scripts/run_benchmarks.sh` orchestrates the full workflow:

1. Stage static inputs to `results/staging/` and build manifests in `results/inputs/`.
2. Generate missing-data simulated VCFs at **5%, 10%, 20%, and 40%** missingness.
3. Run **resource-usage** benchmark on no-missing simulated VCFs (`piawka` vs `pixy`).
4. Run **accuracy/missingness** benchmark on simulated missing-data VCFs.
5. Run **real-data data-gain** benchmark (`piawka` vs `pixy` on the same windows/groups).
6. Run **multiallelic** benchmark (`piawka --mult`) on simulated and real data.
7. Aggregate all run outputs into analysis-ready tables with `scripts/analyze_benchmarks.py`.
8. Build manuscript figures from tables with `scripts/make_figures.py`.

## Repository layout

- `env/`: conda environment specification.
- `scripts/`: executable workflow (`run_benchmarks.sh`), tool wrappers, aggregation, and plotting.
- `results/`: staged inputs, raw benchmark outputs, logs, and derived tables.
  - `results/resource_usage/`
  - `results/accuracy/`
  - `results/data_gain/`
  - `results/multiallelic_check/`
  - `results/tables/`
- `figures/`: generated main and supplementary figures.
- `manuscript/`: paper stub and full methods text describing the executed processing.
- `tools/piawka/`: local `piawka` binary used by wrappers.

## Main outputs

- **Tables:** `results/tables/*.tsv` (resource timing, missingness/theta summaries, data-gain tables, multiallelic pair tables).
- **Figures:** `figures/main/*.png`, `figures/supplementary/*.png`.
- **Per-run raw outputs:** tool-specific outputs and timing logs under `results/resource_usage/`, `results/accuracy/`, `results/data_gain/`, and `results/multiallelic_check/`.

## Reproduce from staged data

```bash
bash scripts/bootstrap_tools.sh
conda run -n piawka-paper-py bash scripts/run_benchmarks.sh
conda run -n piawka-paper-py python scripts/make_figures.py
```

If raw benchmark outputs already exist, rerunning aggregation/plotting only is usually enough:

```bash
conda run -n piawka-paper-py python scripts/analyze_benchmarks.py
conda run -n piawka-paper-py python scripts/make_figures.py
```
