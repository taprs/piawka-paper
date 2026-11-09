#!/usr/bin/env python3
"""Leave-one-South-out jackknife (7 South individuals -> 7 folds; scripts/run_a3_jackknife.sh) of the gene-wise
piawka/pixy Fst percentiles shown in Figure 3. For each fold the percentile is recomputed over the genes defined
in both tools in that fold (100 * (1 - (rank-1)/(N-1)), as in analyze_northsouth_fst_ranks.py). Jackknife
SE = sqrt((n-1)/n * sum_i (theta_i - mean(theta_i))^2), n = 7; 95% CI = full-sample percentile +/- t(0.975, n-1) * SE,
clipped to [0, 100].

Also per set: the delta column (piawka - pixy percentile), delta_percentile (percentile of that delta among all genes
defined in both tools) and per-gene site counts (piawka_snps, pixy_snps).

Output: results/tables/northsouth_jackknife_ci_{gwas_snp,go2000028}_genes.tsv (full-sample columns of the
        northsouth_fst_ranks_* table plus *_jk_se, *_ci_lo, *_ci_hi per tool)
        results/tables/northsouth_jackknife_ci_all_genes.tsv (same CIs for every gene in northsouth_fst_all_genes.tsv,
        matched on chr+end, plus ci_overlap = piawka/pixy CIs overlap)
"""
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = pathlib.Path(__file__).resolve().parents[1]
JK = ROOT / "results" / "a3_jackknife"
TABLES = ROOT / "results" / "tables"
SETS = ("gwas_snp", "go2000028")


def pct(s):
    r = s.rank(ascending=False, method="min")
    return 100 * (1 - (r - 1) / (len(s) - 1))


def load(d):
    pw = pd.read_csv(d / "piawka.bed", sep="\t").rename(columns={"#chr": "chr"})
    pw = pw[(pw["stat"] == "fst") & (pw["locus"] != "-")][["chr", "start", "end", "value"]].rename(columns={"value": "piawka"})
    px = pd.read_csv(d / "pixy" / "pixy_fst.txt", sep="\t")
    px = px.rename(columns={"chromosome": "chr", "window_pos_2": "end", "avg_hudson_fst": "pixy"})
    px["start"] = px["window_pos_1"] - 1
    m = pw.merge(px[["chr", "start", "end", "pixy"]], on=["chr", "start", "end"]).dropna(subset=["piawka", "pixy"])
    m = m.drop_duplicates(["chr", "start", "end"])
    m["interval"] = m["chr"] + ":" + m["start"].astype(str) + "-" + m["end"].astype(str)
    m = m.set_index("interval")
    return pd.DataFrame({"piawka": pct(m["piawka"]), "pixy": pct(m["pixy"])})


folds = {p.name: load(p) for p in sorted(JK.glob("jk*")) if (p / "done").exists()}
n = len(folds)
assert n == 7, f"expected 7 finished folds, found {n}"
print(n, "folds;", {k: len(v) for k, v in folds.items()})
tcrit = t.ppf(0.975, n - 1)

# Per-gene site counts of the full-sample (all individuals) run: piawka segr2 = segregating sites (multiallelic included); pixy no_snps = biallelic SNPs used. Keyed by interval.
STAGING = ROOT / "results" / "staging" / "real"
_pw = pd.read_csv(STAGING / "northsouth_piawka.bed", sep="\t").rename(columns={"#chr": "chr"})
_pw = _pw[(_pw["locus"] != "-") & (_pw["stat"] == "segr2")]
_pw["interval"] = _pw["chr"] + ":" + _pw["start"].astype(str) + "-" + _pw["end"].astype(str)
SITES_PW = _pw.pivot_table(index="interval", columns="stat", values="value", aggfunc="first").rename(
    columns={"segr2": "piawka_snps"})
_px = pd.read_csv(STAGING / "northsouth_pixy_fst.txt", sep="\t")
_px["interval"] = _px["chromosome"] + ":" + (_px["window_pos_1"] - 1).astype(str) + "-" + _px["window_pos_2"].astype(str)
SITES_PX = _px.drop_duplicates("interval").set_index("interval")["no_snps"].rename("pixy_snps")

# Percentile of each gene's delta (piawka percentile - pixy percentile) among all genes defined in both tools.
_all = pd.read_csv(TABLES / "northsouth_fst_all_genes.tsv", sep="\t")
DELTA_PCT = dict(zip(_all["chr"] + ":" + _all["end"].astype(str), pct(_all["delta_pct"]).round(3)))

for s in SETS:
    base = pd.read_csv(TABLES / f"northsouth_fst_ranks_{s}_genes.tsv", sep="\t", keep_default_na=False)
    for tool in ("piawka", "pixy"):
        reps = pd.DataFrame({k: v[tool].reindex(base["interval"]).to_numpy() for k, v in folds.items()})
        if reps.isna().any().any():
            print(f"  {s}/{tool}: {int(reps.isna().any(axis=1).sum())} genes missing in >=1 fold (CI left NA)")
        se = np.sqrt((n - 1) / n * ((reps.sub(reps.mean(axis=1), axis=0)) ** 2).sum(axis=1))
        est = base[f"{tool}_percentile"].astype(float)
        base[f"{tool}_jk_se"] = se.to_numpy()
        base[f"{tool}_ci_lo"] = (est - tcrit * se).clip(lower=0).to_numpy()
        base[f"{tool}_ci_hi"] = (est + tcrit * se).clip(upper=100).to_numpy()
    base = base.rename(columns={"piawka_minus_pixy_percentile": "delta"})
    base["delta_percentile"] = [DELTA_PCT.get(i.split(":")[0] + ":" + i.rsplit("-", 1)[1]) for i in base["interval"]]
    base = base.join(SITES_PW, on="interval").join(SITES_PX, on="interval")
    cols = ["gene_lyrata", "interval", "gene_thaliana", "common_name"]
    for tool in ("piawka", "pixy"):
        cols += [f"{tool}_{c}" for c in ("fst", "percentile", "jk_se", "ci_lo", "ci_hi", "snps")]
    base = base[cols + ["delta", "delta_percentile"]]
    base.to_csv(TABLES / f"northsouth_jackknife_ci_{s}_genes.tsv", sep="\t", index=False, float_format="%.4g")
    w = base[["piawka_jk_se", "pixy_jk_se"]]
    print(f"{s}: {len(base)} genes; median jackknife SE piawka {w.piawka_jk_se.median():.2f}, pixy {w.pixy_jk_se.median():.2f}; "
          f"max {w.max().max():.1f}")

# All genes: northsouth_fst_all_genes.tsv has no start column, so match on chr+end (checked unique below).
allg = pd.read_csv(TABLES / "northsouth_fst_all_genes.tsv", sep="\t")
allg["key"] = allg["chr"] + ":" + allg["end"].astype(str)
assert allg["key"].is_unique
fold_all = {}
for k, v in folds.items():
    v = v.copy()
    v["key"] = [i.split(":")[0] + ":" + i.rsplit("-", 1)[1] for i in v.index]
    assert v["key"].is_unique
    fold_all[k] = v.set_index("key")
for tool in ("piawka", "pixy"):
    reps = pd.DataFrame({k: v[tool].reindex(allg["key"]).to_numpy() for k, v in fold_all.items()})
    print(f"  all/{tool}: {int(reps.isna().any(axis=1).sum())} genes missing in >=1 fold (CI left NA)")
    se = np.sqrt((n - 1) / n * ((reps.sub(reps.mean(axis=1, skipna=False), axis=0)) ** 2).sum(axis=1, skipna=False))
    est = allg[f"{tool}_pct"]
    allg[f"{tool}_jk_se"] = se.to_numpy()
    allg[f"{tool}_ci_lo"] = (est - tcrit * se).clip(lower=0)
    allg[f"{tool}_ci_hi"] = (est + tcrit * se).clip(upper=100)
allg["ci_overlap"] = ~((allg.piawka_ci_lo > allg.pixy_ci_hi) | (allg.pixy_ci_lo > allg.piawka_ci_hi))
allg.drop(columns="key").to_csv(TABLES / "northsouth_jackknife_ci_all_genes.tsv", sep="\t", index=False, float_format="%.4g")
ok = allg.piawka_ci_lo.notna() & allg.pixy_ci_lo.notna()
print(f"all: {len(allg)} genes, {int(ok.sum())} with CIs; non-overlapping CIs: {int((~allg.ci_overlap[ok]).sum())}")
