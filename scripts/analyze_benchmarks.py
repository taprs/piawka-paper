#!/usr/bin/env python3
import csv
import pathlib
import re
import statistics
from collections import defaultdict


ROOT = pathlib.Path(__file__).resolve().parents[1]
BENCH = ROOT / "results"
TABLES = BENCH / "tables"
TABLES.mkdir(parents=True, exist_ok=True)


def read_tsv(path: pathlib.Path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def parse_time_files(base_dir: pathlib.Path):
    rows = []
    for tf in sorted(base_dir.rglob("*.time.tsv")):
        run_id = tf.stem.replace(".time", "")
        m = re.search(r"_t(\d+)$", run_id)
        threads = int(m.group(1)) if m else 1
        for row in read_tsv(tf):
            rows.append(
                {
                    "sample_id": tf.parent.name,
                    "run_id": run_id,
                    "threads": threads,
                    "tool": row.get("tool", ""),
                    "elapsed_sec": row.get("elapsed_sec", ""),
                    "cpu_total_sec": row.get("cpu_total_sec", ""),
                    "max_rss_kb": row.get("max_rss_kb", ""),
                    "time_file": str(tf),
                }
            )
    return rows


def read_piawka_pi(path: pathlib.Path):
    pi = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8 or parts[6] != "pi":
                continue
            chrom, start0, end = parts[0], int(parts[1]), int(parts[2])
            key = (chrom, start0, end, parts[4])
            pi[key] = float(parts[7])
    return pi


def read_piawka_lines(path: pathlib.Path):
    lines = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8 or parts[6] != "lines":
                continue
            chrom, start0, end = parts[0], int(parts[1]), int(parts[2])
            key = (chrom, start0, end, parts[4])
            lines[key] = float(parts[7])
    return lines


def read_piawka_stat(path: pathlib.Path, stat: str):
    vals = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8 or parts[6] != stat:
                continue
            chrom, start0, end = parts[0], int(parts[1]), int(parts[2])
            key = (chrom, start0, end, parts[4])
            vals[key] = float(parts[7])
    return vals


def read_pixy_pi(path: pathlib.Path):
    out = {}
    rows = read_tsv(path)
    for r in rows:
        key = (r["chromosome"], int(r["window_pos_1"]), int(r["window_pos_2"]), r["pop"])
        out[key] = {
            "pi": float(r["avg_pi"]),
            "no_sites": float(r.get("no_sites", 0) or 0),
        }
    return out


def compare_pi(piawka_pi, pixy_pi):
    keys = sorted(set(piawka_pi).intersection(pixy_pi))
    if not keys:
        return None
    abs_diff = []
    rel_diff = []
    for k in keys:
        p = piawka_pi[k]
        x = pixy_pi[k]["pi"]
        d = abs(p - x)
        abs_diff.append(d)
        rel_diff.append(d / x if x else 0.0)
    return {
        "n_overlap": len(keys),
        "mean_abs_diff": statistics.fmean(abs_diff),
        "max_abs_diff": max(abs_diff),
        "mean_rel_diff": statistics.fmean(rel_diff),
        "max_rel_diff": max(rel_diff),
    }


def build_multiallelic_pairs(bial_pi, mult_pi, bial_lines=None, mult_lines=None, pixy_pi=None):
    rows = []
    for chrom, s, e, pop in sorted(set(bial_pi).intersection(mult_pi)):
        bial_val = bial_pi[(chrom, s, e, pop)]
        mult_val = mult_pi[(chrom, s, e, pop)]
        expected = bial_val + (4.0 / 3.0) * (bial_val ** 3)
        row = {
            "chromosome": chrom,
            "window_pos_1": s,
            "window_pos_2": e,
            "group": pop,
            "bial_pi": bial_val,
            "mult_pi": mult_val,
            "expected_mult_pi": expected,
            "abs_diff_bial": abs(mult_val - bial_val),
            "abs_diff_expected": abs(mult_val - expected),
        }
        if bial_lines is not None and (chrom, s, e, pop) in bial_lines:
            row["bial_sites"] = bial_lines[(chrom, s, e, pop)]
        if mult_lines is not None and (chrom, s, e, pop) in mult_lines:
            row["mult_sites"] = mult_lines[(chrom, s, e, pop)]
        if "bial_sites" in row and "mult_sites" in row:
            row["sites_gained_mult_vs_bial"] = row["mult_sites"] - row["bial_sites"]
        if pixy_pi is not None and (chrom, s, e, pop) in pixy_pi:
            row["pixy_pi"] = pixy_pi[(chrom, s, e, pop)]["pi"]
            row["pixy_sites"] = pixy_pi[(chrom, s, e, pop)]["no_sites"]
        rows.append(row)
    return rows


def main():
    # 1) Resource usage timing table
    timing_rows = parse_time_files(BENCH / "resource_usage")
    write_tsv(
        TABLES / "resource_usage_timing.tsv",
        timing_rows,
        ["sample_id", "run_id", "threads", "tool", "elapsed_sec", "cpu_total_sec", "max_rss_kb", "time_file"],
    )
    by_tool_threads = defaultdict(list)
    for r in timing_rows:
        by_tool_threads[(r["tool"], int(r.get("threads", 1)))].append(r)
    timing_thread_rows = []
    for (tool, threads), rows in sorted(by_tool_threads.items(), key=lambda x: (x[0][1], x[0][0])):
        elapsed = [float(x["elapsed_sec"]) for x in rows if x.get("elapsed_sec") not in ("", "NA")]
        rss = [float(x["max_rss_kb"]) for x in rows if x.get("max_rss_kb") not in ("", "NA")]
        timing_thread_rows.append(
            {
                "tool": tool,
                "threads": threads,
                "n_runs": len(rows),
                "elapsed_sec_mean": statistics.fmean(elapsed) if elapsed else "",
                "elapsed_sec_median": statistics.median(elapsed) if elapsed else "",
                "max_rss_kb_mean": statistics.fmean(rss) if rss else "",
                "max_rss_kb_median": statistics.median(rss) if rss else "",
            }
        )
    write_tsv(
        TABLES / "resource_usage_by_threads.tsv",
        timing_thread_rows,
        ["tool", "threads", "n_runs", "elapsed_sec_mean", "elapsed_sec_median", "max_rss_kb_mean", "max_rss_kb_median"],
    )

    # 2) Simulated multiallelic pi comparison: pixy (resource usage) vs piawka --mult
    sim_pairs = []
    for piawka_file in sorted((BENCH / "multiallelic_check" / "simulated_mult").glob("*/piawka_mult.piawka.bed")):
        sample = piawka_file.parent.name
        pixy_candidates = list((BENCH / "resource_usage" / sample / "pixy.pixy").glob("*_pi.txt"))
        if not pixy_candidates:
            continue
        pixy_dict = read_pixy_pi(pixy_candidates[0])
        bial_file = BENCH / "resource_usage" / sample / "piawka.piawka.bed"
        if not bial_file.exists():
            continue
        bial_pi = read_piawka_pi(bial_file)
        bial_lines = read_piawka_lines(bial_file)
        mult_pi = read_piawka_pi(piawka_file)
        mult_lines = read_piawka_lines(piawka_file)
        pairs = build_multiallelic_pairs(bial_pi, mult_pi, bial_lines, mult_lines, pixy_dict)
        for row in pairs:
            sim_pairs.append({"sample_id": sample, **row})
    write_tsv(
        TABLES / "multiallelic_pi_simulated_pairs.tsv",
        sim_pairs,
        [
            "sample_id",
            "chromosome",
            "window_pos_1",
            "window_pos_2",
            "group",
            "bial_pi",
            "mult_pi",
            "bial_sites",
            "mult_sites",
            "sites_gained_mult_vs_bial",
            "pixy_pi",
            "pixy_sites",
            "expected_mult_pi",
            "abs_diff_bial",
            "abs_diff_expected",
        ],
    )
    by_sample = defaultdict(list)
    for row in sim_pairs:
        by_sample[row["sample_id"]].append(row)
    sim_rows = []
    for sample, rows in sorted(by_sample.items()):
        cmp = compare_pi(
            {(r["chromosome"], int(r["window_pos_1"]), int(r["window_pos_2"]), r["group"]): float(r["mult_pi"]) for r in rows},
            {(r["chromosome"], int(r["window_pos_1"]), int(r["window_pos_2"]), r["group"]): {"pi": float(r["bial_pi"])} for r in rows},
        )
        if not cmp:
            continue
        sim_rows.append({"sample_id": sample, **cmp})
        write_tsv(
            TABLES / "multiallelic_pi_simulated.tsv",
            sim_rows,
            ["sample_id", "n_overlap", "mean_abs_diff", "max_abs_diff", "mean_rel_diff", "max_rel_diff"],
        )

    # 3) Real-data per-window pi and site-gain comparisons
    real_piawka_file = BENCH / "data_gain" / "real_lyrata" / "piawka.piawka.bed"
    real_pixy_pi_files = list((BENCH / "data_gain" / "real_lyrata" / "pixy.pixy").glob("*_pi.txt"))
    if real_piawka_file.exists() and real_pixy_pi_files:
        piawka_pi = read_piawka_pi(real_piawka_file)
        piawka_lines = read_piawka_lines(real_piawka_file)
        pixy_pi = read_pixy_pi(real_pixy_pi_files[0])

        # Window-level aligned pi values for line plots.
        pi_window_rows = []
        for k in sorted(set(piawka_pi).intersection(pixy_pi)):
            chrom, s, e, pop = k
            pi_window_rows.append(
                {
                    "chromosome": chrom,
                    "window_pos_1": s,
                    "window_pos_2": e,
                    "group": pop,
                    "piawka_pi": piawka_pi[k],
                    "pixy_pi": pixy_pi[k]["pi"],
                }
            )
        write_tsv(
            TABLES / "data_gain_pi_by_window.tsv",
            pi_window_rows,
            ["chromosome", "window_pos_1", "window_pos_2", "group", "piawka_pi", "pixy_pi"],
        )

        cmp = compare_pi(piawka_pi, pixy_pi)
        if cmp:
            write_tsv(
                TABLES / "multiallelic_pi_real.tsv",
                [cmp],
                ["n_overlap", "mean_abs_diff", "max_abs_diff", "mean_rel_diff", "max_rel_diff"],
            )
        # Real-data multiallelic pairs: pixy (data gain) vs piawka --mult.
        real_mult_file = BENCH / "multiallelic_check" / "real_mult" / "piawka_mult.piawka.bed"
        if real_mult_file.exists():
            real_mult_pairs = build_multiallelic_pairs(
                piawka_pi,
                read_piawka_pi(real_mult_file),
                piawka_lines,
                read_piawka_lines(real_mult_file),
                pixy_pi,
            )
            real_summary = compare_pi(
                {(r["chromosome"], int(r["window_pos_1"]), int(r["window_pos_2"]), r["group"]): float(r["mult_pi"]) for r in real_mult_pairs},
                {(r["chromosome"], int(r["window_pos_1"]), int(r["window_pos_2"]), r["group"]): {"pi": float(r["bial_pi"])} for r in real_mult_pairs},
            )
            if real_summary:
                write_tsv(
                    TABLES / "multiallelic_pi_real.tsv",
                    [real_summary],
                    ["n_overlap", "mean_abs_diff", "max_abs_diff", "mean_rel_diff", "max_rel_diff"],
                )
            write_tsv(
                TABLES / "multiallelic_pi_real_pairs.tsv",
                real_mult_pairs,
                [
                    "chromosome",
                    "window_pos_1",
                    "window_pos_2",
                    "group",
                    "bial_pi",
                    "mult_pi",
                    "bial_sites",
                    "mult_sites",
                    "sites_gained_mult_vs_bial",
                    "pixy_pi",
                    "pixy_sites",
                    "expected_mult_pi",
                    "abs_diff_bial",
                    "abs_diff_expected",
                ],
            )

        # Per-group site-count comparison (piawka lines vs pixy no_sites)
        site_rows = []
        gain_rows = []
        for k in sorted(set(piawka_lines).intersection(pixy_pi)):
            chrom, s, e, pop = k
            pixy_pi_val = pixy_pi[k]["pi"]
            piawka_pi_val = piawka_pi.get(k)
            site_gain = piawka_lines[k] - pixy_pi[k]["no_sites"]
            site_rows.append(
                {
                    "chromosome": chrom,
                    "window_pos_1": s,
                    "window_pos_2": e,
                    "group": pop,
                    "piawka_lines": piawka_lines[k],
                    "pixy_no_sites": pixy_pi[k]["no_sites"],
                    "sites_gained_piawka_vs_pixy": site_gain,
                }
            )
            if piawka_pi_val is not None:
                site_gain_pct = ""
                if pixy_pi[k]["no_sites"]:
                    site_gain_pct = (site_gain / pixy_pi[k]["no_sites"]) * 100.0
                gain_rows.append(
                    {
                        "chromosome": chrom,
                        "window_pos_1": s,
                        "window_pos_2": e,
                        "group": pop,
                        "pixy_pi": pixy_pi_val,
                        "piawka_pi": piawka_pi_val,
                        "pi_gain_piawka_vs_pixy": piawka_pi_val - pixy_pi_val,
                        "pixy_sites": pixy_pi[k]["no_sites"],
                        "piawka_sites": piawka_lines[k],
                        "sites_gained_piawka_vs_pixy": site_gain,
                        "sites_gained_pct_vs_pixy": site_gain_pct,
                    }
                )
        write_tsv(
            TABLES / "data_gain_sites_per_window.tsv",
            site_rows,
            [
                "chromosome",
                "window_pos_1",
                "window_pos_2",
                "group",
                "piawka_lines",
                "pixy_no_sites",
                "sites_gained_piawka_vs_pixy",
            ],
        )
        write_tsv(
            TABLES / "data_gain_pi_vs_site_gain.tsv",
            gain_rows,
            [
                "chromosome",
                "window_pos_1",
                "window_pos_2",
                "group",
                "pixy_pi",
                "piawka_pi",
                "pi_gain_piawka_vs_pixy",
                "pixy_sites",
                "piawka_sites",
                "sites_gained_piawka_vs_pixy",
                "sites_gained_pct_vs_pixy",
            ],
        )

        # 100-fold leave-one-window-out fold means and SD by group/tool
        by_group_tool = defaultdict(lambda: {"piawka": {}, "pixy": {}})
        for k, val in piawka_pi.items():
            chrom, s, e, pop = k
            by_group_tool[pop]["piawka"][(chrom, s, e)] = val
        for k, val in pixy_pi.items():
            chrom, s, e, pop = k
            by_group_tool[pop]["pixy"][(chrom, s, e)] = val["pi"]

        cv_fold_rows = []
        cv_sd_rows = []
        for pop, d in by_group_tool.items():
            windows = sorted(set(d["piawka"]).intersection(d["pixy"]))
            if len(windows) < 2:
                continue
            folds_p = []
            folds_x = []
            for i, w in enumerate(windows, start=1):
                pvals = [d["piawka"][ww] for ww in windows if ww != w]
                xvals = [d["pixy"][ww] for ww in windows if ww != w]
                pmean = statistics.fmean(pvals)
                xmean = statistics.fmean(xvals)
                folds_p.append(pmean)
                folds_x.append(xmean)
                cv_fold_rows.append(
                    {"group": pop, "tool": "piawka", "fold_id": i, "fold_mean_pi": pmean}
                )
                cv_fold_rows.append(
                    {"group": pop, "tool": "pixy", "fold_id": i, "fold_mean_pi": xmean}
                )
            cv_sd_rows.append(
                {"group": pop, "tool": "piawka", "n_folds": len(folds_p), "mean_of_folds": statistics.fmean(folds_p), "sd_of_folds": statistics.stdev(folds_p)}
            )
            cv_sd_rows.append(
                {"group": pop, "tool": "pixy", "n_folds": len(folds_x), "mean_of_folds": statistics.fmean(folds_x), "sd_of_folds": statistics.stdev(folds_x)}
            )
        write_tsv(
            TABLES / "data_gain_crossval_folds.tsv",
            cv_fold_rows,
            ["group", "tool", "fold_id", "fold_mean_pi"],
        )
        write_tsv(
            TABLES / "data_gain_crossval_sd.tsv",
            cv_sd_rows,
            ["group", "tool", "n_folds", "mean_of_folds", "sd_of_folds"],
        )

    # 4) Accuracy benchmark: theta_w comparison and pi-minus-theta proxies under missingness.
    baseline_rows = []
    for sample_dir in sorted((BENCH / "resource_usage").glob("*")):
        pfile = sample_dir / "piawka.piawka.bed"
        x_pi_file = sample_dir / "pixy.pixy" / "pixy_pi.txt"
        x_tw_file = sample_dir / "pixy.pixy" / "pixy_watterson_theta.txt"
        if not (pfile.exists() and x_pi_file.exists() and x_tw_file.exists()):
            continue
        p_pi = read_piawka_stat(pfile, "pi")
        p_tw = read_piawka_stat(pfile, "theta_w")
        x_pi = read_pixy_pi(x_pi_file)
        x_tw_rows = read_tsv(x_tw_file)
        x_tw = {}
        for r in x_tw_rows:
            k = (r["chromosome"], int(r["window_pos_1"]), int(r["window_pos_2"]), r["pop"])
            x_tw[k] = float(r["avg_watterson_theta"])
        keys = sorted(set(p_pi).intersection(p_tw).intersection(x_pi).intersection(x_tw))
        if not keys:
            continue
        k = keys[0]
        baseline_rows.append(
            {
                "sample_id": sample_dir.name,
                "source_class": "",
                "piawka_pi_base": p_pi[k],
                "pixy_pi_base": x_pi[k]["pi"],
                "piawka_theta_w_base": p_tw[k],
                "pixy_theta_w_base": x_tw[k],
                "piawka_pi_minus_theta_w_base": p_pi[k] - p_tw[k],
                "pixy_pi_minus_theta_w_base": x_pi[k]["pi"] - x_tw[k],
                "piawka_theta_low_base": "",
                "piawka_pi_minus_theta_low_base": "",
            }
        )
    baseline_by_sample = {r["sample_id"]: r for r in baseline_rows}

    source_class_by_sample = {}
    manifest_nomiss = BENCH / "inputs" / "sim_nomissing_manifest.tsv"
    if manifest_nomiss.exists():
        for r in read_tsv(manifest_nomiss):
            source_class_by_sample[r["sample_id"]] = r["source_class"]
    for r in baseline_rows:
        r["source_class"] = source_class_by_sample.get(r["sample_id"], "")
        # No-missing theta_low baseline is read from accuracy/<sample>/r000 if available.
        p0 = BENCH / "accuracy" / r["sample_id"] / "r000" / "piawka.piawka.bed"
        if p0.exists():
            p0_pi = read_piawka_stat(p0, "pi")
            p0_tl = read_piawka_stat(p0, "theta_low")
            keys0 = sorted(set(p0_pi).intersection(p0_tl))
            if keys0:
                k0 = keys0[0]
                r["piawka_theta_low_base"] = p0_tl[k0]
                r["piawka_pi_minus_theta_low_base"] = p0_pi[k0] - p0_tl[k0]

    acc_rows = []
    # Add explicit no-missing baseline rows from resource benchmark.
    for base in baseline_rows:
        acc_rows.append(
            {
                "sample_id": base["sample_id"],
                "source_class": base.get("source_class", ""),
                "missing_rate": 0.0,
                "piawka_pi": base["piawka_pi_base"],
                "pixy_pi": base["pixy_pi_base"],
                "piawka_theta_w": base["piawka_theta_w_base"],
                "pixy_theta_w": base["pixy_theta_w_base"],
                "piawka_theta_low": base["piawka_theta_low_base"],
                "piawka_pi_minus_theta_w": base["piawka_pi_minus_theta_w_base"],
                "pixy_pi_minus_theta_w": base["pixy_pi_minus_theta_w_base"],
                "piawka_pi_minus_theta_low": base["piawka_pi_minus_theta_low_base"],
                "piawka_theta_w_rel_to_nomiss": 1.0,
                "pixy_theta_w_rel_to_nomiss": 1.0,
                "piawka_pi_minus_theta_w_rel_to_nomiss": 1.0,
                "pixy_pi_minus_theta_w_rel_to_nomiss": 1.0,
                "piawka_pi_minus_theta_low_rel_to_nomiss": 1.0 if base["piawka_pi_minus_theta_low_base"] != "" else "",
            }
        )
    for sample_dir in sorted((BENCH / "accuracy").glob("*")):
        if not sample_dir.is_dir():
            continue
        sample_id = sample_dir.name
        for rate_dir in sorted(sample_dir.glob("r*")):
            pfile = rate_dir / "piawka.piawka.bed"
            x_pi_file = rate_dir / "pixy.pixy" / "pixy_pi.txt"
            x_tw_file = rate_dir / "pixy.pixy" / "pixy_watterson_theta.txt"
            if not (pfile.exists() and x_pi_file.exists() and x_tw_file.exists()):
                continue
            # rate directory tags are r005/r010/r020/r040 => 0.05/0.10/0.20/0.40
            rate = int(rate_dir.name[1:]) / 100.0
            p_pi = read_piawka_stat(pfile, "pi")
            p_tw = read_piawka_stat(pfile, "theta_w")
            p_tl = read_piawka_stat(pfile, "theta_low")
            x_pi = read_pixy_pi(x_pi_file)
            x_tw = {}
            for r in read_tsv(x_tw_file):
                k = (r["chromosome"], int(r["window_pos_1"]), int(r["window_pos_2"]), r["pop"])
                x_tw[k] = float(r["avg_watterson_theta"])
            keys = sorted(set(p_pi).intersection(p_tw).intersection(p_tl).intersection(x_pi).intersection(x_tw))
            if not keys:
                continue
            k = keys[0]
            row = {
                "sample_id": sample_id,
                "source_class": source_class_by_sample.get(sample_id, ""),
                "missing_rate": rate,
                "piawka_pi": p_pi[k],
                "pixy_pi": x_pi[k]["pi"],
                "piawka_theta_w": p_tw[k],
                "pixy_theta_w": x_tw[k],
                "piawka_theta_low": p_tl[k],
            }
            row["piawka_pi_minus_theta_w"] = row["piawka_pi"] - row["piawka_theta_w"]
            row["pixy_pi_minus_theta_w"] = row["pixy_pi"] - row["pixy_theta_w"]
            row["piawka_pi_minus_theta_low"] = row["piawka_pi"] - row["piawka_theta_low"]
            base = baseline_by_sample.get(sample_id)
            if base:
                row["piawka_theta_w_rel_to_nomiss"] = row["piawka_theta_w"] / base["piawka_theta_w_base"] if base["piawka_theta_w_base"] else ""
                row["pixy_theta_w_rel_to_nomiss"] = row["pixy_theta_w"] / base["pixy_theta_w_base"] if base["pixy_theta_w_base"] else ""
                row["piawka_pi_minus_theta_w_rel_to_nomiss"] = row["piawka_pi_minus_theta_w"] / base["piawka_pi_minus_theta_w_base"] if base["piawka_pi_minus_theta_w_base"] else ""
                row["pixy_pi_minus_theta_w_rel_to_nomiss"] = row["pixy_pi_minus_theta_w"] / base["pixy_pi_minus_theta_w_base"] if base["pixy_pi_minus_theta_w_base"] else ""
                row["piawka_pi_minus_theta_low_rel_to_nomiss"] = row["piawka_pi_minus_theta_low"] / base["piawka_pi_minus_theta_low_base"] if base["piawka_pi_minus_theta_low_base"] else ""
            acc_rows.append(row)

    write_tsv(
        TABLES / "accuracy_theta_proxy_per_run.tsv",
        acc_rows,
        [
            "sample_id",
            "source_class",
            "missing_rate",
            "piawka_pi",
            "pixy_pi",
            "piawka_theta_w",
            "pixy_theta_w",
            "piawka_theta_low",
            "piawka_pi_minus_theta_w",
            "pixy_pi_minus_theta_w",
            "piawka_pi_minus_theta_low",
            "piawka_theta_w_rel_to_nomiss",
            "pixy_theta_w_rel_to_nomiss",
            "piawka_pi_minus_theta_w_rel_to_nomiss",
            "pixy_pi_minus_theta_w_rel_to_nomiss",
            "piawka_pi_minus_theta_low_rel_to_nomiss",
        ],
    )

    by_rate = defaultdict(list)
    for r in acc_rows:
        by_rate[r["missing_rate"]].append(r)
    acc_summary = []
    for rate in sorted(by_rate):
        rows = by_rate[rate]
        def avg(field):
            vals = [x[field] for x in rows if field in x and x[field] != ""]
            return statistics.fmean(vals) if vals else ""
        acc_summary.append(
            {
                "missing_rate": rate,
                "n": len(rows),
                "piawka_theta_w_mean": avg("piawka_theta_w"),
                "pixy_theta_w_mean": avg("pixy_theta_w"),
                "piawka_pi_minus_theta_w_mean": avg("piawka_pi_minus_theta_w"),
                "pixy_pi_minus_theta_w_mean": avg("pixy_pi_minus_theta_w"),
                "piawka_pi_minus_theta_low_mean": avg("piawka_pi_minus_theta_low"),
                "piawka_theta_w_rel_to_nomiss_mean": avg("piawka_theta_w_rel_to_nomiss"),
                "pixy_theta_w_rel_to_nomiss_mean": avg("pixy_theta_w_rel_to_nomiss"),
                "piawka_pi_minus_theta_w_rel_to_nomiss_mean": avg("piawka_pi_minus_theta_w_rel_to_nomiss"),
                "pixy_pi_minus_theta_w_rel_to_nomiss_mean": avg("pixy_pi_minus_theta_w_rel_to_nomiss"),
                "piawka_pi_minus_theta_low_rel_to_nomiss_mean": avg("piawka_pi_minus_theta_low_rel_to_nomiss"),
            }
        )
    write_tsv(
        TABLES / "accuracy_theta_proxy_summary.tsv",
        acc_summary,
        [
            "missing_rate",
            "n",
            "piawka_theta_w_mean",
            "pixy_theta_w_mean",
            "piawka_pi_minus_theta_w_mean",
            "pixy_pi_minus_theta_w_mean",
            "piawka_pi_minus_theta_low_mean",
            "piawka_theta_w_rel_to_nomiss_mean",
            "pixy_theta_w_rel_to_nomiss_mean",
            "piawka_pi_minus_theta_w_rel_to_nomiss_mean",
            "pixy_pi_minus_theta_w_rel_to_nomiss_mean",
            "piawka_pi_minus_theta_low_rel_to_nomiss_mean",
        ],
    )

    # Category-level distributions for sinaplot visualization (onepop vs bneck).
    cat_rows = []
    by_cat = defaultdict(list)
    for r in acc_rows:
        if not r.get("source_class"):
            continue
        by_cat[r["source_class"]].append(r)
    for cat, rows in sorted(by_cat.items()):
        def avg(field):
            vals = [x[field] for x in rows if field in x and x[field] != ""]
            return statistics.fmean(vals) if vals else ""
        cat_rows.append(
            {
                "source_class": cat,
                "n": len(rows),
                "piawka_pi_minus_theta_w_mean": avg("piawka_pi_minus_theta_w"),
                "pixy_pi_minus_theta_w_mean": avg("pixy_pi_minus_theta_w"),
                "piawka_pi_minus_theta_low_mean": avg("piawka_pi_minus_theta_low"),
            }
        )
    write_tsv(
        TABLES / "accuracy_theta_proxy_by_category.tsv",
        cat_rows,
        [
            "source_class",
            "n",
            "piawka_pi_minus_theta_w_mean",
            "pixy_pi_minus_theta_w_mean",
            "piawka_pi_minus_theta_low_mean",
        ],
    )

    # Category-by-missingness summary for grouped sinaplot reporting.
    cat_miss_rows = []
    by_cat_rate = defaultdict(list)
    for r in acc_rows:
        if not r.get("source_class"):
            continue
        by_cat_rate[(r["source_class"], r["missing_rate"])].append(r)
    for (cat, rate), rows in sorted(by_cat_rate.items(), key=lambda x: (x[0][0], x[0][1])):
        def avg(field):
            vals = [x[field] for x in rows if field in x and x[field] != ""]
            return statistics.fmean(vals) if vals else ""
        cat_miss_rows.append(
            {
                "source_class": cat,
                "missing_rate": rate,
                "n": len(rows),
                "piawka_pi_minus_theta_w_mean": avg("piawka_pi_minus_theta_w"),
                "pixy_pi_minus_theta_w_mean": avg("pixy_pi_minus_theta_w"),
                "piawka_pi_minus_theta_low_mean": avg("piawka_pi_minus_theta_low"),
            }
        )
    write_tsv(
        TABLES / "accuracy_theta_proxy_by_category_missingness.tsv",
        cat_miss_rows,
        [
            "source_class",
            "missing_rate",
            "n",
            "piawka_pi_minus_theta_w_mean",
            "pixy_pi_minus_theta_w_mean",
            "piawka_pi_minus_theta_low_mean",
        ],
    )

    print("Wrote benchmark summary tables in results/tables")


if __name__ == "__main__":
    main()
