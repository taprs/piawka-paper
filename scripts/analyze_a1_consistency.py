#!/usr/bin/env python3
"""A1 (review response): piawka vs pixy agreement on the 16 rescaled two-population SLiM replicates,
without and with simulated missing data.
Reads results/a1_twopop16/<Rep>/{none,masked}/{default,biallelic}/{piawka.bed,pixy/*.txt} (scripts/run_a1_twopop.sh).
Writes results/tables/a1_twopop_per_rep.tsv (one row per replicate x missing x variant x stat x population) and
results/tables/a1_twopop_summary.tsv (n, Pearson r, Lin's CCC, mean and max |relative difference|, mean site counts).
"""
import csv
import itertools
import os
import pathlib
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
A1 = ROOT / "results" / "a1_twopop16"
A1P = A1
SUF = ""
TABLES = ROOT / "results" / "tables"
MISS = ("none", "masked")
VARIANTS = ("default", "biallelic")


def tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def piawka(path):
    out, sites = {}, {}
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        stat, pop1, pop2, val = f[6], f[4], f[5], float(f[7])
        if stat == "lines":
            sites[pop1] = val
        else:
            out[(stat, pop1 if pop2 == "." else "p1-p2")] = val
    out[("fst", "p1-p2")] = out.pop(("fst", "p1-p2"))
    return out, max(sites.values())


def pixy(d):
    out, sites = {}, None
    for r in tsv(d / "pixy_pi.txt"):
        out[("pi", r["pop"])] = float(r["avg_pi"]); sites = float(r["no_sites"])
    for r in tsv(d / "pixy_dxy.txt"):
        out[("dxy", "p1-p2")] = float(r["avg_dxy"])
    for r in tsv(d / "pixy_fst.txt"):
        out[("fst", "p1-p2")] = float(r["avg_hudson_fst"])
    for r in tsv(d / "pixy_watterson_theta.txt"):
        out[("theta_w", r["pop"])] = float(r["avg_watterson_theta"])
    for r in tsv(d / "pixy_tajima_d.txt"):
        out[("tajima", r["pop"])] = float(r["tajima_d"])
    return out, sites


def ccc(x, y):
    mx, my = x.mean(), y.mean()
    return 2 * np.cov(x, y, bias=True)[0, 1] / (x.var() + y.var() + (mx - my) ** 2)


rows = []
for rep_dir in sorted(A1.glob("Rep*"), key=lambda p: int(p.name[3:])):
    if not (rep_dir / "done").exists():
        continue
    for miss, variant in itertools.product(MISS, VARIANTS):
        d = rep_dir / miss / variant
        pw, pw_sites = piawka(A1P / rep_dir.name / miss / variant / "piawka.bed")
        px, px_sites = pixy(d / "pixy")
        for (stat, pop), v in pw.items():
            if stat == "fstwc" or (stat, pop) not in px:
                continue
            rows.append(dict(rep=rep_dir.name, missing=miss, variant=variant, stat=stat, pop=pop, piawka=v, pixy=px[(stat, pop)],
                             abs_diff=v - px[(stat, pop)], rel_diff_pct=100 * (v - px[(stat, pop)]) / abs(px[(stat, pop)]),
                             sites_piawka=pw_sites, sites_pixy=px_sites))

cols = ["rep", "missing", "variant", "stat", "pop", "piawka", "pixy", "abs_diff", "rel_diff_pct", "sites_piawka", "sites_pixy"]
with (TABLES / f"a1_twopop_per_rep{SUF}.tsv").open("w") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)

summ = []
for miss, variant in itertools.product(MISS, VARIANTS):
    for stat in ("pi", "dxy", "fst", "theta_w", "tajima"):
        sel = [r for r in rows if r["missing"] == miss and r["variant"] == variant and r["stat"] == stat]
        x = np.array([r["piawka"] for r in sel]); y = np.array([r["pixy"] for r in sel])
        rel = (x - y) / np.abs(y)
        summ.append(dict(missing=miss, variant=variant, stat=stat, n=len(sel), pearson_r=np.corrcoef(x, y)[0, 1], lin_ccc=ccc(x, y),
                         mean_rel_diff_pct=100 * rel.mean(), max_abs_rel_diff_pct=100 * np.abs(rel).max(),
                         max_abs_diff=np.abs(x - y).max(),
                         mean_sites_piawka=np.mean([r["sites_piawka"] for r in sel]), mean_sites_pixy=np.mean([r["sites_pixy"] for r in sel])))
with (TABLES / f"a1_twopop_summary{SUF}.tsv").open("w") as fh:
    w = csv.DictWriter(fh, list(summ[0]), delimiter="\t"); w.writeheader(); w.writerows(summ)
for s in summ:
    print("\t".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in s.values()))
