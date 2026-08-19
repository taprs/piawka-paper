#!/usr/bin/env python3
import csv
import pathlib
import re


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


def build_multiallelic_pairs(bial_pi, mult_pi, bial_lines=None, mult_lines=None, pixy_pi=None, pixy_mult_pi=None):
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
        if pixy_mult_pi is not None and (chrom, s, e, pop) in pixy_mult_pi:
            row["pixymult_pi"] = pixy_mult_pi[(chrom, s, e, pop)]["pi"]
            row["pixymult_sites"] = pixy_mult_pi[(chrom, s, e, pop)]["no_sites"]
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
    # 1b) Window-size comparison at a single process: one 10 Mbp window
    #     (syn_win.bed) vs eight 1.25 Mbp windows (syn_win_threads8.bed).
    win_rows = []
    for sample_dir in sorted((BENCH / "resource_usage_compare_1thread").glob("*")):
        if not sample_dir.is_dir():
            continue
        for tool in ("piawka", "pixy"):
            base_f = sample_dir / f"{tool}_syn_win.time.tsv"
            cmp_f = sample_dir / f"{tool}_syn_win_threads8.time.tsv"
            if not (base_f.exists() and cmp_f.exists()):
                continue
            base = read_tsv(base_f)[0]
            cmp_ = read_tsv(cmp_f)[0]

            def num(row, field):
                v = row.get(field, "")
                return float(v) if v not in ("", "NA") else None

            base_t, cmp_t = num(base, "elapsed_sec"), num(cmp_, "elapsed_sec")
            base_m, cmp_m = num(base, "max_rss_kb"), num(cmp_, "max_rss_kb")
            win_rows.append(
                {
                    "sample_id": sample_dir.name,
                    "tool": tool,
                    "threads": 1,
                    "window_set_baseline": "syn_win.bed",
                    "window_set_compare": "syn_win_threads8.bed",
                    "elapsed_sec_syn_win": base_t if base_t is not None else "",
                    "elapsed_sec_syn_win_threads8": cmp_t if cmp_t is not None else "",
                    "elapsed_ratio_threads8_over_syn_win": (cmp_t / base_t) if base_t else "",
                    "max_rss_kb_syn_win": base_m if base_m is not None else "",
                    "max_rss_kb_syn_win_threads8": cmp_m if cmp_m is not None else "",
                    "max_rss_ratio_threads8_over_syn_win": (cmp_m / base_m) if base_m else "",
                }
            )
    write_tsv(
        TABLES / "resource_usage_1thread_synwin_vs_threads8.tsv",
        win_rows,
        [
            "sample_id",
            "tool",
            "threads",
            "window_set_baseline",
            "window_set_compare",
            "elapsed_sec_syn_win",
            "elapsed_sec_syn_win_threads8",
            "elapsed_ratio_threads8_over_syn_win",
            "max_rss_kb_syn_win",
            "max_rss_kb_syn_win_threads8",
            "max_rss_ratio_threads8_over_syn_win",
        ],
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
    # 3) Real-data per-window pi and site-gain comparisons
    real_piawka_file = BENCH / "data_gain" / "real_lyrata" / "piawka.piawka.bed"
    real_pixy_pi_files = list((BENCH / "data_gain" / "real_lyrata" / "pixy.pixy").glob("*_pi.txt"))
    if real_piawka_file.exists() and real_pixy_pi_files:
        piawka_pi = read_piawka_pi(real_piawka_file)
        piawka_lines = read_piawka_lines(real_piawka_file)
        pixy_pi = read_pixy_pi(real_pixy_pi_files[0])

        # Real-data multiallelic pairs: pixy (data gain) vs piawka --mult, plus pixy --include_multiallelic_snps.
        real_mult_file = BENCH / "multiallelic_check" / "real_mult" / "piawka_mult.piawka.bed"
        real_pixy_mult_pi_files = list(
            (BENCH / "multiallelic_check" / "real_mult" / "pixy_mult.pixy").glob("*_pi.txt")
        )
        pixy_mult_pi = read_pixy_pi(real_pixy_mult_pi_files[0]) if real_pixy_mult_pi_files else None
        if real_mult_file.exists():
            real_mult_pairs = build_multiallelic_pairs(
                piawka_pi,
                read_piawka_pi(real_mult_file),
                piawka_lines,
                read_piawka_lines(real_mult_file),
                pixy_pi,
                pixy_mult_pi,
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
                    "pixymult_pi",
                    "pixymult_sites",
                    "expected_mult_pi",
                    "abs_diff_bial",
                    "abs_diff_expected",
                ],
            )

        # Per-group site-count and pi comparison (piawka lines vs pixy no_sites)
        gain_rows = []
        for k in sorted(set(piawka_lines).intersection(pixy_pi)):
            chrom, s, e, pop = k
            pixy_pi_val = pixy_pi[k]["pi"]
            piawka_pi_val = piawka_pi.get(k)
            site_gain = piawka_lines[k] - pixy_pi[k]["no_sites"]
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

    print("Wrote benchmark summary tables in results/tables")


if __name__ == "__main__":
    main()
