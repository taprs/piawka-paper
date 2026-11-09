#!/usr/bin/env python3
"""Collect the resource-usage runs into one tidy table, results/tables/a2_scaling.tsv (Figure 1, Supplementary Figure S1):
real data (results/a2_real, scripts/run_a2_real_scaling.sh) and the 16 simulated neutral replicates
(results/resource_sim16, scripts/run_resource_sim16.sh). One row per tool x run with series (real_samples, real_groups,
sim_small, sim_windows, sim_threads), n_samples, n_groups, processes, window_mbp, records, elapsed_sec, max_rss_mb.
Also collects the profiling table (results/a2_profile) and prints linear fits and speedups."""
import pathlib
import re

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
rows = []


def add(d, series, tool, **kw):
    t = d / f"{tool}.time.tsv"
    if not t.exists():
        return
    r = pd.read_csv(t, sep="\t").iloc[0]
    rec = (d / "records.txt").read_text().strip() if (d / "records.txt").exists() else kw.pop("records", np.nan)
    rows.append(dict(series=series, tool=tool, records=float(rec), elapsed_sec=r.elapsed_sec,
                     cpu_total_sec=r.cpu_total_sec, max_rss_mb=r.max_rss_kb / 1024, **kw))


REAL_RECORDS = 1074832
for d in sorted((ROOT / "results/a2_real").glob("*_p*")):
    m = re.match(r"(n|k)(\d+)_p(\d+)", d.name)
    if not m or not (d / "done").exists():
        continue
    kind, v, p = m.group(1), int(m.group(2)), int(m.group(3))
    for tool in ("piawka", "pixy"):
        add(d, "real_samples" if kind == "n" else "real_groups", tool,
            n_samples=v if kind == "n" else 910, n_groups=12 if kind == "n" else v,
            processes=p, window_mbp=0.01, records=REAL_RECORDS)
# simulated arm: 16 rescaled neutral one-population replicates (scripts/run_resource_sim16.sh); 100 samples = win1250000_p{1,8}
for d in sorted((ROOT / "results/resource_sim16").glob("*/*_p*")):
    if not (d / "done").exists():
        continue
    rep = int(re.search(r"Rep(\d+)", d.parent.name).group(1))
    m = re.match(r"n(\d+)_p(\d+)", d.name)
    w = re.match(r"win(\d+)_p(\d+)", d.name)
    for tool in ("piawka", "pixy"):
        if m:
            add(d, "sim_small", tool, replicate=rep, n_samples=int(m.group(1)), n_groups=1,
                processes=int(m.group(2)), window_mbp=1.25)
        elif w:
            wm, p = int(w.group(1)) / 1e6, int(w.group(2))
            if p == 1:
                add(d, "sim_windows", tool, replicate=rep, n_samples=100, n_groups=1, processes=p, window_mbp=wm)
            if wm == 1.25:  # eight windows, 100 samples: the process series (Fig 1B) and the 100-sample point of the sample series
                add(d, "sim_threads", tool, replicate=rep, n_samples=100, n_groups=1, processes=p, window_mbp=wm)
                if p in (1, 8):
                    add(d, "sim_small", tool, replicate=rep, n_samples=100, n_groups=1, processes=p, window_mbp=wm)
df = pd.DataFrame(rows)
df.to_csv(ROOT / "results/tables/a2_scaling.tsv", sep="\t", index=False, float_format="%.4g")
print(df.groupby(["series", "tool", "processes"]).size().to_string())
for (s, tool, p), g in df[df.series.isin(["real_samples", "sim_small"])].groupby(["series", "tool", "processes"]):
    m = g.groupby("n_samples")[["elapsed_sec", "max_rss_mb"]].median()
    if len(m) > 2:
        r_t = np.corrcoef(m.index, m.elapsed_sec)[0, 1]
        r_m = np.corrcoef(m.index, m.max_rss_mb)[0, 1]
        print(f"{s} {tool} p{p}: r(time,n)={r_t:.4f} r(rss,n)={r_m:.4f}")

# profiling stage ablation (scripts/run_a2_profile.sh): median of 3 repeats per stage
pf = ROOT / "results/a2_profile/profile.tsv"
if pf.exists():
    p = pd.read_csv(pf, sep="\t")
    p["stage"] = p["stage"].str[0]
    p = p.groupby(["stage", "description"], as_index=False).elapsed_sec.median()
    p.to_csv(ROOT / "results/tables/a2_profile.tsv", sep="\t", index=False, float_format="%.3f")
    print(p.to_string())
