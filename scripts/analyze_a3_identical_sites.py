#!/usr/bin/env python3
"""A3.2 (review response): gene-wise North/South FST, piawka 27e6ae3 vs pixy 2.2.3, default site filtering vs both tools on the
same biallelic SNPs. Reads results/a3_northsouth/<variant>/{piawka.bed,pixy/pixy_fst.txt} (scripts/run_a3_northsouth.sh).
Writes results/tables/northsouth_identical_sites.tsv (one row per variant)."""
import csv, pathlib
import numpy as np
from scipy import stats
ROOT = pathlib.Path(__file__).resolve().parents[1]
A3 = ROOT / "results" / "a3_northsouth"

def load(variant):
    pk = {}
    for r in csv.reader(open(A3 / variant / "piawka.bed"), delimiter="\t"):
        if r[0].startswith("#") or r[6] != "fst" or r[7] in ("NA", "nan", ""): continue
        pk[(r[0], int(r[2]))] = float(r[7])
    px = {}
    for r in csv.DictReader(open(A3 / variant / "pixy" / "pixy_fst.txt"), delimiter="\t"):
        if r["avg_hudson_fst"] in ("NA", ""): continue
        px[(r["chromosome"], int(r["window_pos_2"]))] = float(r["avg_hudson_fst"])
    return pk, px

def pct(x):
    rk = stats.rankdata(-x, method="min"); return 100 * (1 - (rk - 1) / (len(x) - 1))

rows = []
for variant in ("default", "biallelic"):
    pk, px = load(variant)
    keys = sorted(set(pk) & set(px))
    a = np.array([pk[k] for k in keys]); b = np.array([px[k] for k in keys]); d = pct(a) - pct(b)
    rows.append(dict(variant=variant, n_piawka=len(pk), n_pixy=len(px), n_both=len(keys),
                     spearman=stats.spearmanr(a, b)[0], max_abs_fst_diff=np.abs(a - b).max(),
                     median_abs_fst_diff=np.median(np.abs(a - b)),
                     median_abs_delta_pct=np.median(np.abs(d)), share_abs_delta_gt10=np.mean(np.abs(d) > 10),
                     share_abs_delta_gt25=np.mean(np.abs(d) > 25)))
with (ROOT / "results" / "tables" / "northsouth_identical_sites.tsv").open("w") as fh:
    w = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
for r in rows: print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
