#!/usr/bin/env python3
"""Bundle the supplementary tables into one .xlsx workbook, in the order of the analyses, plus a column legend:
  1 resource use (a2_*.tsv)                         2 piawka vs pixy consistency, 16 two-population replicates, with and without
                                                      simulated missing data (scripts/analyze_a1_consistency.py)
  3, 4 North/South gene-wise Fst rank tables (GWAS SNP genes, GO2000028 genes; scripts/analyze_northsouth_fst_ranks.py, annotated with
       leave-one-South-out jackknife 95% CIs, percentile of the piawka-pixy delta and per-gene site counts by scripts/analyze_a3_jackknife.py)
  5 piawka-only sites in the genes where piawka and pixy percentiles differ (scripts/analyze_a3_skew_sites.py); observed ALT allele(s) underlined.

Usage:
  conda run -n piawka-paper python scripts/collect_fst_tables.py [output.xlsx]

Defaults to results/tables/piawka_supplementary_data.xlsx.
"""
import csv
import pathlib
import statistics
import sys

from openpyxl import Workbook
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont

ROOT = pathlib.Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"

SHEETS = [
    ("Supp Table 3 - GWAS SNP genes", TABLES / "northsouth_jackknife_ci_gwas_snp_genes.tsv"),
    ("Supp Table 4 - GO2000028 genes", TABLES / "northsouth_jackknife_ci_go2000028_genes.tsv"),
]

LEGEND = [
    ("piawka_percentile / pixy_percentile", "Gene-wise Fst percentile among all genes defined in both tools (100 = highest Fst)"),
    ("delta", "piawka percentile minus pixy percentile"),
    ("delta_percentile", "Percentile of delta among all genes defined in both tools (100 = piawka percentile most above pixy's)"),
    ("*_jk_se, *_ci_lo, *_ci_hi", "Leave-one-South-out jackknife (7 folds) SE and 95% CI (t, 6 df) of the percentile, clipped to [0, 100]"),
    ("piawka_snps", "piawka: segregating sites in the window, multiallelic included (segr2)"),
    ("pixy_snps", "pixy: biallelic SNPs used in the window (no_snps)"),
]

# Supp Table 5: sites piawka uses but pixy cannot, in the genes with non-overlapping jackknife CIs
# (scripts/analyze_a3_skew_sites.py). The ALT allele(s) observed in North/South are underlined in the alt column
# (explained in the caption, not the legend sheet); the helper column alt_observed is not written to the sheet.
SKEW_SHEET = ("Supp Table 5 - skew sites", TABLES / "northsouth_skew_sites_piawka_only.tsv")

# Supp Table 2: every value compared in the consistency check (one row per replicate x missing x variant x stat x population)
CONSISTENCY_SHEET = ("Supp Table 2 - consistency", TABLES / "a1_twopop_per_rep.tsv")
LEGEND += [
    ("rep", "Simulated two-population replicate (Demotwopop_scaled10_Rep0-15)"),
    ("missing", "none: simulated VCF as is; masked: each genotype set to missing with a per-site probability drawn from U(0.2, 0.6), seeded"),
    ("variant", "default: each tool's own site filtering; biallelic: both tools on the same VCF restricted to biallelic SNPs (bcftools view -m2 -M2 -v snps)"),
    ("stat / pop", "pi, theta_w and tajima per population (p1, p2); dxy and fst between populations (p1-p2); fst is Hudson's"),
    ("piawka / pixy", "Value reported by each tool for the 10 Mbp window"),
    ("abs_diff / rel_diff_pct", "piawka minus pixy, absolute and as % of the pixy value (relative differences of Tajima's D are inflated where D is near 0)"),
    ("sites_piawka / sites_pixy", "Sites used by each tool"),
]

# Supp Table 1: resource use of piawka and pixy (a2_scaling.tsv, a2_profile.tsv, a2_window_memory_fit.tsv), runs of the
# same condition summarised as mean and SD (simulated series: 16 replicates; real data: one run per condition)
RESOURCE_SHEET = "Supp Table 1 - resource use"
RESOURCE_SERIES = [
    ("sim_windows", "Simulated: window size (10 Mbp analysed as windows of 0.625-10 Mbp, 1 process)"),
    ("sim_threads", "Simulated: number of processes (eight 1.25 Mbp windows)"),
    ("sim_small", "Simulated: number of samples (10-75 samples, 1.25 Mbp windows)"),
    ("real_samples", "Real A. lyrata: number of samples"),
    ("real_groups", "Real A. lyrata: number of populations (910 samples)"),
]
LEGEND += [
    ("n_runs", "Runs of the condition: 16 simulated replicates, or a single run on the real data"),
    ("elapsed_sec / cpu_total_sec / max_rss_mb", "Wall-clock time, total CPU time (all processes) and peak resident memory (mean and SD over runs)"),
    ("window_mbp", "Window size in Mbp (0.01 = the 10 kb windows of the real-data runs)"),
    ("profiling rows", "piawka stage ablation on 200 kb of the real data (tabix only, + field splitting, + genotype scan, then piawka statistics)"),
    ("window-size fit rows", "Linear fit of peak memory on window size (simulated, 1 process): slope in MB per Mbp with 95% CI, R2"),
]

DEFAULT_OUTPUT = TABLES / "piawka_supplementary_data.xlsx"


# decimals per numeric column of Supp Table 2: fixed-point, no scientific notation (site counts are integers)
CONSISTENCY_DECIMALS = {"piawka": 6, "pixy": 6, "abs_diff": 10, "rel_diff_pct": 8}


def _txt(v, decimals):
    """Numeric TSV field as fixed-point text, so Excel does not reformat it."""
    try:
        x = float(v)
    except ValueError:
        return v
    return str(int(x)) if decimals is None else f"{x:.{decimals}f}"


def main():
    output = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUTPUT

    wb = Workbook()
    wb.remove(wb.active)

    for sheet_name, tsv_path in SHEETS:
        if not tsv_path.exists():
            raise FileNotFoundError(f"missing input table: {tsv_path}")
        ws = wb.create_sheet(title=sheet_name)
        with tsv_path.open(newline="") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                ws.append(row)
        for col_cells in ws.columns:
            width = max(len(str(c.value)) for c in col_cells if c.value is not None)
            ws.column_dimensions[col_cells[0].column_letter].width = min(width + 2, 40)

    sheet_name, tsv_path = SKEW_SHEET
    if not tsv_path.exists():
        raise FileNotFoundError(f"missing input table: {tsv_path}")
    ws = wb.create_sheet(title=sheet_name)
    with tsv_path.open(newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header, body = rows[0], rows[1:]
    i_alt, i_obs = header.index("alt"), header.index("alt_observed")
    ws.append(header[:i_obs] + header[i_obs + 1:])
    for row in body:
        observed = set(row[i_obs].split(","))
        parts = []
        for k, allele in enumerate(row[i_alt].split(",")):
            if k:
                parts.append(",")
            parts.append(TextBlock(InlineFont(u="single"), allele) if allele in observed else allele)
        values = list(row)  # everything as text, so signs and trailing zeros survive
        values[i_alt] = CellRichText(*parts)
        ws.append(values[:i_obs] + values[i_obs + 1:])
    for col_cells in ws.columns:
        width = max(len(str(c.value)) for c in col_cells if c.value is not None)
        ws.column_dimensions[col_cells[0].column_letter].width = min(width + 2, 40)

    sheet_name, tsv_path = CONSISTENCY_SHEET
    if not tsv_path.exists():
        raise FileNotFoundError(f"missing input table: {tsv_path}")
    ws = wb.create_sheet(title=sheet_name)
    with tsv_path.open(newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
        decimals = [CONSISTENCY_DECIMALS.get(c) for c in rows[0]]
        ws.append(rows[0])
        for row in rows[1:]:
            ws.append([_txt(v, d) if c in CONSISTENCY_DECIMALS or c.startswith("sites_") else v
                       for v, d, c in zip(row, decimals, rows[0])])
    for col_cells in ws.columns:
        width = max(len(str(c.value)) for c in col_cells if c.value is not None)
        ws.column_dimensions[col_cells[0].column_letter].width = min(width + 2, 24)

    ws = wb.create_sheet(title=RESOURCE_SHEET)
    with (TABLES / "a2_scaling.tsv").open(newline="") as fh:
        runs = list(csv.DictReader(fh, delimiter="\t"))
    ws.append(["benchmark", "tool", "n_samples", "n_populations", "processes", "window_mbp", "records", "n_runs",
               "elapsed_sec_mean", "elapsed_sec_sd", "cpu_total_sec_mean", "cpu_total_sec_sd", "max_rss_mb_mean", "max_rss_mb_sd"])
    for series, label in RESOURCE_SERIES:
        groups = {}
        for r in runs:
            if r["series"] == series:
                groups.setdefault((r["n_samples"], r["n_groups"], r["processes"], r["window_mbp"], r["tool"]), []).append(r)
        for key in sorted(groups, key=lambda k: (float(k[0]), float(k[1]), float(k[2]), float(k[3]), k[4])):
            g = groups[key]
            n_samples, n_groups, processes, window, tool = key
            row = [label, tool, int(float(n_samples)), int(float(n_groups)), int(float(processes)), float(window),
                   int(float(g[0]["records"])), len(g)]
            for col in ("elapsed_sec", "cpu_total_sec", "max_rss_mb"):
                v = [float(r[col]) for r in g]
                row += [round(statistics.mean(v), 3), round(statistics.stdev(v), 3) if len(v) > 1 else None]
            ws.append(row)
    ws.append([])
    ws.append(["piawka profiling (200 kb, real data)", "stage", "elapsed_sec"])
    with (TABLES / "a2_profile.tsv").open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ws.append(["", r["description"], float(r["elapsed_sec"])])
    ws.append([])
    ws.append(["window-size fit (simulated, 1 process)", "tool", "n", "slope_mb_per_mbp", "slope_ci_lo", "slope_ci_hi", "r2"])
    with (TABLES / "a2_window_memory_fit.tsv").open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ws.append(["", r["tool"], int(r["n"])] + [round(float(r[c]), 4) for c in ("slope_mb_per_mbp", "slope_ci_lo", "slope_ci_hi", "r2")])
    ws.column_dimensions["A"].width = 80
    for col in "BCDEFGHIJKLMN":
        ws.column_dimensions[col].width = 18

    wb._sheets.sort(key=lambda w: (w.title == "Column legend", w.title))
    ws = wb.create_sheet(title="Column legend")
    ws.append(["column", "description"])
    for row in LEGEND:
        ws.append(list(row))
    ws.column_dimensions["A"].width = 40
    ws.column_dimensions["B"].width = 110

    output.parent.mkdir(parents=True, exist_ok=True)
    wb.save(output)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
