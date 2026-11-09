#!/usr/bin/env python3
import csv
import pathlib
import re


ROOT = pathlib.Path(__file__).resolve().parents[1]
import os
BENCH = ROOT / "results" / os.environ.get("RESULTS_SUFFIX", "")  # RESULTS_SUFFIX: rerun output tree (see run_benchmarks.sh)
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
        # pixy >= 2.2 reports 1-based window starts; piawka BED starts are 0-based (pixy 2.0.0.beta14 reported them 0-based)
        key = (r["chromosome"], int(r["window_pos_1"]) - 1, int(r["window_pos_2"]), r["pop"])
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
    # 1) Simulated multiallelic pi comparison (16 neutral replicates): piawka --mult (run_benchmarks.sh) vs the biallelic
    #    piawka and pixy runs over the same 10 Mbp window (run_resource_sim16.sh, results/resource_sim16/<Rep>/win10000000_p1)
    sim_pairs = []
    for piawka_file in sorted((BENCH / "multiallelic_check" / "simulated_mult").glob("*/piawka_mult.piawka.bed")):
        sample = piawka_file.parent.name
        base = ROOT / "results" / "resource_sim16" / sample / "win10000000_p1"
        pixy_candidates = list((base / "pixy.pixy").glob("*_pi.txt"))
        if not pixy_candidates:
            continue
        pixy_dict = read_pixy_pi(pixy_candidates[0])
        bial_file = base / "piawka.piawka.bed"
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
    # 2) Real-data per-window pi and site-gain comparisons
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
