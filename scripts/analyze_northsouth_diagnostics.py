#!/usr/bin/env python3
"""Genome-wide diagnostics of the North/South FST piawka-vs-pixy comparison (existing outputs only).

Inputs : results/staging/real/{northsouth_piawka.bed,northsouth_pixy_fst.txt}, results/tables/northsouth_fst_ranks_*.tsv
Outputs: results/tables/northsouth_fst_diagnostics.tsv (one test/statistic per row)
         results/tables/northsouth_fst_all_genes.tsv   (per-gene values, percentiles, pixy SNP count)
Percentiles follow analyze_northsouth_fst_ranks.py (100 = most differentiated, over genes defined in both tools).
"""
import csv, pathlib
import numpy as np
from scipy import stats
ROOT = pathlib.Path(__file__).resolve().parents[1]
ST = ROOT / "results" / "staging" / "real"; TB = ROOT / "results" / "tables"
pk = {}
for r in csv.reader(open(ST / "northsouth_piawka.bed"), delimiter="\t"):
    if r[0].startswith("#") or r[6] != "fst" or r[7] in ("NA", "nan", ""): continue
    pk[(r[0], int(r[2]))] = (r[3], float(r[7]), float(r[8]), float(r[9]), int(r[1]))
px = {}
for r in csv.DictReader(open(ST / "northsouth_pixy_fst.txt"), delimiter="\t"):
    if r["avg_hudson_fst"] in ("NA", ""): continue
    px[(r["chromosome"], int(r["window_pos_2"]))] = (float(r["avg_hudson_fst"]), int(r["no_snps"]))
keys = sorted(set(pk) & set(px))
print("piawka", len(pk), "pixy", len(px), "matched", len(keys))
a = np.array([pk[k][1] for k in keys]); b = np.array([px[k][0] for k in keys]); snps = np.array([px[k][1] for k in keys])
def pct(x):
    rk = stats.rankdata(-x, method="min"); n = len(x); return 100 * (1 - (rk - 1) / (n - 1))
pa, pb = pct(a), pct(b); d = pa - pb
rows = []
def add(t, v, p=""): rows.append((t, v, p)); print(t, v, p)
add("n_genes_matched", len(keys))
rho = stats.spearmanr(a, b); add("spearman_fst_all_genes", round(rho[0], 4), rho[1])
add("delta_pct_median_all_genes", round(float(np.median(d)), 3))
add("delta_pct_iqr_all_genes", "%.2f..%.2f" % tuple(np.percentile(d, [25, 75])))
add("delta_pct_p5_p95_all_genes", "%.2f..%.2f" % tuple(np.percentile(d, [5, 95])))
add("share_abs_delta_gt10_all_genes", round(float(np.mean(np.abs(d) > 10)), 4))
add("share_abs_delta_gt25_all_genes", round(float(np.mean(np.abs(d) > 25)), 4))
r2 = stats.spearmanr(snps, np.abs(d)); add("spearman_abs_delta_vs_pixy_no_snps", round(r2[0], 4), r2[1])
lo = snps <= np.percentile(snps, 25); add("median_abs_delta_lowest_quartile_pixy_snps", round(float(np.median(np.abs(d[lo]))), 3))
add("median_abs_delta_other_quartiles", round(float(np.median(np.abs(d[~lo]))), 3))
den = np.array([pk[k][3] for k in keys]); r3 = stats.spearmanr(den, np.abs(d)); add("spearman_abs_delta_vs_piawka_denominator", round(r3[0], 4), r3[1])
for name in ("gwas_snp", "go2000028"):
    t = list(csv.DictReader(open(TB / f"northsouth_fst_ranks_{name}_genes.tsv"), delimiter="\t"))
    dd = np.array([float(x["piawka_percentile"]) - float(x["pixy_percentile"]) for x in t])
    w = stats.wilcoxon(dd); add(f"{name}_n", len(t)); add(f"{name}_wilcoxon_paired_delta", round(float(w.statistic), 2), w.pvalue)
    for tool in ("piawka", "pixy"):
        v = np.array([float(x[f"{tool}_percentile"]) for x in t])
        k = int((v >= 80).sum()); bt = stats.binomtest(k, len(v), 0.2, alternative="greater")
        mw = stats.mannwhitneyu(v, pa if tool == "piawka" else pb, alternative="greater")
        add(f"{name}_{tool}_top20_count", f"{k}/{len(v)}", bt.pvalue)
        add(f"{name}_{tool}_median_pct", round(float(np.median(v)), 2), mw.pvalue)
with open(TB / "northsouth_fst_diagnostics.tsv", "w") as f:
    f.write("quantity\tvalue\tp_value\n"); [f.write("\t".join(map(str, r)) + "\n") for r in rows]
with open(TB / "northsouth_fst_all_genes.tsv", "w") as f:
    f.write("chr\tend\tlocus\tpiawka_fst\tpiawka_den\tpixy_fst\tpixy_no_snps\tpiawka_pct\tpixy_pct\tdelta_pct\n")
    for k, i in zip(keys, range(len(keys))):
        f.write(f"{k[0]}\t{k[1]}\t{pk[k][0]}\t{a[i]}\t{pk[k][3]}\t{b[i]}\t{snps[i]}\t{pa[i]:.3f}\t{pb[i]:.3f}\t{d[i]:.3f}\n")
