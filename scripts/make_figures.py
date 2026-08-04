#!/usr/bin/env python3
import pathlib
import random
import subprocess
import sys

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = pathlib.Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
FIG_MAIN = ROOT / "figures" / "main"
FIG_SUPP = ROOT / "figures" / "supplementary"
FIG_MAIN.mkdir(parents=True, exist_ok=True)
FIG_SUPP.mkdir(parents=True, exist_ok=True)


def group_color_map(groups):
    groups = list(groups)
    n = len(groups)
    if n <= 20:
        cmap = plt.get_cmap("tab20", n)
        colors = [cmap(i) for i in range(n)]
    else:
        cmap = plt.get_cmap("gist_ncar", n)
        colors = [cmap(i) for i in range(n)]
    return {g: colors[i] for i, g in enumerate(groups)}


def main() -> int:
    timing = TABLES / "resource_usage_timing.tsv"
    sim_pairs = TABLES / "multiallelic_pi_simulated_pairs.tsv"
    real_pairs = TABLES / "multiallelic_pi_real_pairs.tsv"
    gain_sites = TABLES / "data_gain_sites_per_window.tsv"
    gain_pi_vs_site_gain = TABLES / "data_gain_pi_vs_site_gain.tsv"
    gain_pi_by_window = TABLES / "data_gain_pi_by_window.tsv"
    acc_proxy = TABLES / "accuracy_theta_proxy_summary.tsv"
    acc_proxy_per_run = TABLES / "accuracy_theta_proxy_per_run.tsv"
    resource_by_threads = TABLES / "resource_usage_by_threads.tsv"
    resource_1thread_compare = TABLES / "resource_usage_1thread_synwin_vs_threads8.tsv"
    resource_1thread_groups_compare = TABLES / "resource_usage_1thread_groups10_vs_full.tsv"
    for req in (
        timing,
        sim_pairs,
        real_pairs,
        gain_sites,
        gain_pi_vs_site_gain,
        gain_pi_by_window,
        acc_proxy,
        acc_proxy_per_run,
        resource_by_threads,
        resource_1thread_compare,
        resource_1thread_groups_compare,
    ):
        if not req.exists():
            raise FileNotFoundError(f"Missing table: {req}")

    # Main figure 1: resource usage scatter (wall-time vs RSS per run).
    tdf = pd.read_csv(timing, sep="\t")
    if tdf.empty:
        raise ValueError("resource_usage_timing.tsv is empty")
    tdf["elapsed_sec"] = pd.to_numeric(tdf["elapsed_sec"], errors="coerce")
    tdf["max_rss_mb"] = pd.to_numeric(tdf["max_rss_kb"], errors="coerce") / 1024.0
    tdf["threads"] = pd.to_numeric(tdf["threads"], errors="coerce")
    tdf = tdf[(tdf["elapsed_sec"] > 0) & (tdf["max_rss_mb"] > 0)]
    fig, ax = plt.subplots(figsize=(7, 5))
    tool_colors = {"piawka": "tab:blue", "pixy": "tab:orange"}
    thread_markers = {1: "o", 2: "s", 4: "^", 8: "D"}

    for tool, sub in tdf.groupby("tool"):
        # Connect the same dataset across thread counts.
        for _, ssub in sub.groupby("sample_id"):
            ssub = ssub.sort_values("threads")
            if len(ssub) > 1:
                ax.plot(
                    ssub["elapsed_sec"],
                    ssub["max_rss_mb"],
                    color=tool_colors.get(tool, "gray"),
                    linewidth=0.7,
                    alpha=0.25,
                )
        for th, th_sub in sub.groupby("threads"):
            marker = thread_markers.get(int(th), "o")
            ax.scatter(
                th_sub["elapsed_sec"],
                th_sub["max_rss_mb"],
                alpha=0.75,
                s=24,
                marker=marker,
                color=tool_colors.get(tool, "gray"),
            )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)")
    ax.set_ylabel("Max RSS (MB)")
    ax.set_title("Resource usage per run")
    tool_handles = [
        Line2D([0], [0], marker="o", linestyle="", color=tool_colors.get(t, "gray"), label=t)
        for t in sorted(tdf["tool"].dropna().unique())
    ]
    thread_handles = [
        Line2D([0], [0], marker=m, linestyle="", color="black", label=f"{int(t)} thread")
        for t, m in sorted(thread_markers.items())
        if t in set(int(x) for x in tdf["threads"].dropna().unique())
    ]
    legend1 = ax.legend(handles=tool_handles, title="Tool", loc="upper left")
    ax.add_artist(legend1)
    ax.legend(handles=thread_handles, title="Threads", loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG_MAIN / "figure1_resource_time_vs_rss.png", dpi=300)
    plt.close(fig)

    # Supplementary: piawka resource usage across thread counts.
    rbt = pd.read_csv(resource_by_threads, sep="\t")
    rbt["threads"] = pd.to_numeric(rbt["threads"], errors="coerce")
    rbt["elapsed_sec_mean"] = pd.to_numeric(rbt["elapsed_sec_mean"], errors="coerce")
    rbt["max_rss_mb_mean"] = pd.to_numeric(rbt["max_rss_kb_mean"], errors="coerce") / 1024.0
    pr = rbt[rbt["tool"] == "piawka"].sort_values("threads")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(pr["threads"], pr["elapsed_sec_mean"], marker="o", linewidth=1.8)
    axes[0].set_xlabel("Threads")
    axes[0].set_ylabel("Mean elapsed time (s)")
    axes[0].set_title("piawka runtime vs threads")
    axes[1].plot(pr["threads"], pr["max_rss_mb_mean"], marker="o", linewidth=1.8, color="tab:red")
    axes[1].set_xlabel("Threads")
    axes[1].set_ylabel("Mean max RSS (MB)")
    axes[1].set_title("piawka memory vs threads")
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "supplementary_piawka_resource_by_threads.png", dpi=300)
    plt.close(fig)

    # Supplementary: 1-thread resource comparison across window files.
    rc = pd.read_csv(resource_1thread_compare, sep="\t")
    for c in ["elapsed_sec_syn_win", "elapsed_sec_syn_win_threads8", "max_rss_kb_syn_win", "max_rss_kb_syn_win_threads8"]:
        rc[c] = pd.to_numeric(rc[c], errors="coerce")
    comp = (
        rc.groupby("tool", as_index=False)[
            ["elapsed_sec_syn_win", "elapsed_sec_syn_win_threads8", "max_rss_kb_syn_win", "max_rss_kb_syn_win_threads8"]
        ]
        .mean()
        .sort_values("tool")
    )
    tools = list(comp["tool"])
    x = range(len(tools))
    w = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar([i - w / 2 for i in x], comp["elapsed_sec_syn_win"], width=w, label="syn_win.bed")
    axes[0].bar([i + w / 2 for i in x], comp["elapsed_sec_syn_win_threads8"], width=w, label="syn_win_threads8.bed")
    axes[0].set_xticks(list(x))
    axes[0].set_xticklabels(tools)
    axes[0].set_ylabel("Mean elapsed time (s)")
    axes[0].set_title("1-thread runtime by window file")
    axes[0].legend()
    axes[1].bar([i - w / 2 for i in x], comp["max_rss_kb_syn_win"] / 1024.0, width=w, label="syn_win.bed")
    axes[1].bar([i + w / 2 for i in x], comp["max_rss_kb_syn_win_threads8"] / 1024.0, width=w, label="syn_win_threads8.bed")
    axes[1].set_xticks(list(x))
    axes[1].set_xticklabels(tools)
    axes[1].set_ylabel("Mean max RSS (MB)")
    axes[1].set_title("1-thread memory by window file")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "supplementary_resource_1thread_window_compare.png", dpi=300)
    plt.close(fig)

    # Supplementary: 1-thread resource comparison across group files.
    rg = pd.read_csv(resource_1thread_groups_compare, sep="\t")
    for c in ["elapsed_sec_full_groups", "elapsed_sec_groups10", "max_rss_kb_full_groups", "max_rss_kb_groups10"]:
        rg[c] = pd.to_numeric(rg[c], errors="coerce")
    comp = (
        rg.groupby("tool", as_index=False)[
            ["elapsed_sec_full_groups", "elapsed_sec_groups10", "max_rss_kb_full_groups", "max_rss_kb_groups10"]
        ]
        .mean()
        .sort_values("tool")
    )
    tools = list(comp["tool"])
    x = range(len(tools))
    w = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar([i - w / 2 for i in x], comp["elapsed_sec_full_groups"], width=w, label="syn_groups.tsv")
    axes[0].bar([i + w / 2 for i in x], comp["elapsed_sec_groups10"], width=w, label="syn_groups_10.tsv")
    axes[0].set_xticks(list(x))
    axes[0].set_xticklabels(tools)
    axes[0].set_ylabel("Mean elapsed time (s)")
    axes[0].set_title("1-thread runtime by groups file")
    axes[0].legend()
    axes[1].bar([i - w / 2 for i in x], comp["max_rss_kb_full_groups"] / 1024.0, width=w, label="syn_groups.tsv")
    axes[1].bar([i + w / 2 for i in x], comp["max_rss_kb_groups10"] / 1024.0, width=w, label="syn_groups_10.tsv")
    axes[1].set_xticks(list(x))
    axes[1].set_xticklabels(tools)
    axes[1].set_ylabel("Mean max RSS (MB)")
    axes[1].set_title("1-thread memory by groups file")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "supplementary_resource_1thread_groups_compare.png", dpi=300)
    plt.close(fig)

    # Main figure 2: simulated biallelic vs multiallelic pi.
    sim_df = pd.read_csv(sim_pairs, sep="\t")
    sim_df["bial_pi"] = pd.to_numeric(sim_df["bial_pi"], errors="coerce")
    sim_df["mult_pi"] = pd.to_numeric(sim_df["mult_pi"], errors="coerce")
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(sim_df["bial_pi"], sim_df["mult_pi"], alpha=0.6, s=14)
    x_max = max(float(sim_df["bial_pi"].max()), float(sim_df["mult_pi"].max()), 1e-6)
    x = pd.Series([0.0, x_max * 1.02])
    y = x
    ax.plot(x, y, color="black", linewidth=1.5, linestyle="--", label=r"$\pi_{mult}=\pi_{bial}$")
    ax.set_xlabel("Biallelic pi (pixy)")
    ax.set_ylabel("Multiallelic pi (piawka --mult)")
    ax.set_title("Simulated multiallelic benchmark")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_MAIN / "figure2_multiallelic_simulated_expectation.png", dpi=300)
    plt.close(fig)

    # Supplementary: simulated pixy vs piawka (both biallelic).
    sim_df["pixy_pi"] = pd.to_numeric(sim_df["pixy_pi"], errors="coerce")
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(sim_df["pixy_pi"], sim_df["bial_pi"], alpha=0.65, s=18)
    x_max = max(float(sim_df["pixy_pi"].max()), float(sim_df["bial_pi"].max()), 1e-6)
    ax.plot([0.0, x_max * 1.02], [0.0, x_max * 1.02], color="black", linewidth=1.5, linestyle="--", label=r"$\pi_{piawka}=\pi_{pixy}$")
    ax.set_xlabel("Biallelic pi (pixy)")
    ax.set_ylabel("Biallelic pi (piawka)")
    ax.set_title("Simulated biallelic comparison")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "supplementary_simulated_pixy_vs_piawka_bial.png", dpi=300)
    plt.close(fig)

    # Main figure 3: real-data piawka biallelic vs piawka multiallelic pi.
    real_df = pd.read_csv(real_pairs, sep="\t")
    real_df["bial_pi"] = pd.to_numeric(real_df["bial_pi"], errors="coerce")
    real_df["mult_pi"] = pd.to_numeric(real_df["mult_pi"], errors="coerce")
    real_groups = sorted(real_df["group"].dropna().unique())
    real_colors = group_color_map(real_groups)
    fig, ax = plt.subplots(figsize=(6, 6))
    for grp in real_groups:
        sub = real_df[real_df["group"] == grp]
        ax.scatter(sub["bial_pi"], sub["mult_pi"], alpha=0.65, s=18, label=grp, color=real_colors[grp])
    x_max = max(float(real_df["bial_pi"].max()), float(real_df["mult_pi"].max()), 1e-6)
    x = pd.Series([0.0, x_max * 1.02])
    y = x
    ax.plot(x, y, color="black", linewidth=1.5, linestyle="--", label=r"$\pi_{mult}=\pi_{bial}$")
    ax.set_xlabel("Biallelic pi (piawka)")
    ax.set_ylabel("Multiallelic pi (piawka --mult)")
    ax.set_title("Real-data multiallelic benchmark")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_MAIN / "figure3_multiallelic_real_expectation.png", dpi=300)
    plt.close(fig)

    # Supplementary: site-count distributions by group.
    sites_df = pd.read_csv(gain_sites, sep="\t")
    sites_df["piawka_lines"] = pd.to_numeric(sites_df["piawka_lines"], errors="coerce")
    sites_df["pixy_no_sites"] = pd.to_numeric(sites_df["pixy_no_sites"], errors="coerce")
    site_groups = sorted(sites_df["group"].dropna().unique())
    site_colors = group_color_map(site_groups)
    fig, ax = plt.subplots(figsize=(8, 5))
    for grp in site_groups:
        sub = sites_df[sites_df["group"] == grp]
        ax.scatter(sub["pixy_no_sites"], sub["piawka_lines"], alpha=0.35, s=10, label=grp, color=site_colors[grp])
    xy_max = float(max(sites_df["pixy_no_sites"].max(), sites_df["piawka_lines"].max(), 1))
    ax.plot([0, xy_max], [0, xy_max], linestyle="--", linewidth=1, color="black")
    ax.set_xlabel("pixy sites")
    ax.set_ylabel("piawka lines")
    ax.set_title("Data gain by group: sites per window")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "supplementary_data_gain_sites.png", dpi=300)
    plt.close(fig)

    # Supplementary: pi gain vs site gain (piawka vs pixy).
    gain_df = pd.read_csv(gain_pi_vs_site_gain, sep="\t")
    gain_df["pi_gain_piawka_vs_pixy"] = pd.to_numeric(gain_df["pi_gain_piawka_vs_pixy"], errors="coerce")
    gain_df["sites_gained_pct_vs_pixy"] = pd.to_numeric(gain_df["sites_gained_pct_vs_pixy"], errors="coerce")
    gain_groups = sorted(gain_df["group"].dropna().unique())
    gain_colors = group_color_map(gain_groups)
    fig, ax = plt.subplots(figsize=(8, 5))
    for grp in gain_groups:
        sub = gain_df[gain_df["group"] == grp]
        ax.scatter(
            sub["sites_gained_pct_vs_pixy"],
            sub["pi_gain_piawka_vs_pixy"],
            alpha=0.5,
            s=14,
            label=grp,
            color=gain_colors[grp],
        )
    ax.axhline(0, color="black", linewidth=1, linestyle="--")
    ax.set_xlabel("Sites gained by piawka vs pixy (%)")
    ax.set_ylabel("pi gain (piawka - pixy)")
    ax.set_title("Data gain benchmark: pi gain vs % sites gain")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "supplementary_data_gain_pi_gain_vs_sites_gain.png", dpi=300)
    plt.close(fig)

    # Supplementary: theta_watterson and proxy stability under missing data.
    ap = pd.read_csv(acc_proxy, sep="\t")
    ap["missing_pct"] = pd.to_numeric(ap["missing_rate"], errors="coerce") * 100.0
    ap = ap.sort_values("missing_pct")
    for c in [
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
    ]:
        ap[c] = pd.to_numeric(ap[c], errors="coerce")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(ap["missing_pct"], ap["piawka_theta_w_mean"], marker="o", label="piawka theta_w")
    axes[0].plot(ap["missing_pct"], ap["pixy_theta_w_mean"], marker="o", label="pixy theta_watterson")
    axes[0].set_xlabel("Introduced missing data (%)")
    axes[0].set_ylabel("Mean theta_watterson")
    axes[0].set_title("Theta_watterson comparison")
    axes[0].legend()

    axes[1].plot(ap["missing_pct"], ap["piawka_pi_minus_theta_w_rel_to_nomiss_mean"], marker="o", label="piawka (pi - theta_w), rel")
    axes[1].plot(ap["missing_pct"], ap["pixy_pi_minus_theta_w_rel_to_nomiss_mean"], marker="o", label="pixy (pi - theta_w), rel")
    axes[1].plot(ap["missing_pct"], ap["piawka_pi_minus_theta_low_rel_to_nomiss_mean"], marker="o", label="piawka (pi - theta_low), rel")
    axes[1].axhline(1.0, color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Introduced missing data (%)")
    axes[1].set_ylabel("Relative mean (stability)")
    axes[1].set_title("Selection/demography proxy stability")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(FIG_SUPP / "supplementary_accuracy_theta_proxy_stability.png", dpi=300)
    plt.close(fig)

    # Supplementary: category-level proxy distributions (sinaplot style).
    apr = pd.read_csv(acc_proxy_per_run, sep="\t")
    apr["missing_rate"] = pd.to_numeric(apr["missing_rate"], errors="coerce")
    for c in ["piawka_pi_minus_theta_w", "pixy_pi_minus_theta_w", "piawka_pi_minus_theta_low"]:
        apr[c] = pd.to_numeric(apr[c], errors="coerce")
    cats = sorted([c for c in apr["source_class"].dropna().unique() if c])
    miss_levels = sorted(apr["missing_rate"].dropna().unique())
    miss_colors = {m: plt.get_cmap("viridis", len(miss_levels))(i) for i, m in enumerate(miss_levels)}
    proxies = [
        ("piawka_pi_minus_theta_w", "piawka: pi - theta_w"),
        ("pixy_pi_minus_theta_w", "pixy: pi - theta_w"),
        ("piawka_pi_minus_theta_low", "piawka: pi - theta_low"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), sharey=False)
    rng = random.Random(42)
    for ax, (field, title) in zip(axes, proxies):
        ticks, labels = [], []
        n_m = len(miss_levels)
        for i, cat in enumerate(cats):
            for j, mr in enumerate(miss_levels):
                sub = apr[(apr["source_class"] == cat) & (apr["missing_rate"] == mr)][field].dropna()
                if len(sub) == 0:
                    continue
                subgroup_center = i + (j - (n_m - 1) / 2.0) * 0.16
                xs = [subgroup_center + rng.uniform(-0.04, 0.04) for _ in range(len(sub))]
                ax.scatter(xs, sub, alpha=0.5, s=12, color=miss_colors[mr])
                mean = float(sub.mean())
                sd = float(sub.std(ddof=1)) if len(sub) > 1 else 0.0
                ax.errorbar([subgroup_center], [mean], yerr=[sd], fmt="o", color="black", capsize=3, markersize=3)
            ticks.append(i)
            labels.append(cat)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_title(title)
        ax.set_xlabel("Simulated data class")
    axes[0].set_ylabel("Proxy value")
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=miss_colors[m], markersize=6, label=f"{int(round(m * 100))}%")
        for m in miss_levels
    ]
    axes[-1].legend(handles=legend_handles, title="Missingness", loc="best")
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "supplementary_accuracy_proxy_category_sinaplot.png", dpi=300)
    plt.close(fig)

    # Supplementary: real-data pixy vs piawka (both biallelic) with same layout.
    pi_df = pd.read_csv(gain_pi_by_window, sep="\t")
    pi_df["piawka_pi"] = pd.to_numeric(pi_df["piawka_pi"], errors="coerce")
    pi_df["pixy_pi"] = pd.to_numeric(pi_df["pixy_pi"], errors="coerce")
    pi_groups = sorted(pi_df["group"].dropna().unique())
    pi_colors = group_color_map(pi_groups)
    fig, ax = plt.subplots(figsize=(6, 6))
    for grp in pi_groups:
        sub = pi_df[pi_df["group"] == grp]
        ax.scatter(sub["pixy_pi"], sub["piawka_pi"], alpha=0.65, s=18, label=grp, color=pi_colors[grp])
    x_max = max(float(pi_df["pixy_pi"].max()), float(pi_df["piawka_pi"].max()), 1e-6)
    ax.plot([0.0, x_max * 1.02], [0.0, x_max * 1.02], color="black", linewidth=1.5, linestyle="--", label=r"$\pi_{piawka}=\pi_{pixy}$")
    ax.set_xlabel("Biallelic pi (pixy)")
    ax.set_ylabel("Biallelic pi (piawka)")
    ax.set_title("Real-data biallelic comparison")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "supplementary_real_pixy_vs_piawka_bial.png", dpi=300)
    plt.close(fig)

    # Supplementary: pixy/piawka pi along chromosome as lines.
    pi_df["window_pos_1"] = pd.to_numeric(pi_df["window_pos_1"], errors="coerce")
    groups = sorted(pi_df["group"].dropna().unique())
    fig, axes = plt.subplots(len(groups), 1, figsize=(10, 2.6 * max(1, len(groups))), sharex=True)
    if len(groups) == 1:
        axes = [axes]
    for ax, grp in zip(axes, groups):
        sub = pi_df[pi_df["group"] == grp].sort_values(["chromosome", "window_pos_1"])
        ax.plot(sub["window_pos_1"], sub["piawka_pi"], label="piawka", linewidth=1.5)
        ax.plot(sub["window_pos_1"], sub["pixy_pi"], label="pixy", linewidth=1.5)
        ax.set_ylabel(f"{grp}\npi")
        ax.legend(loc="upper right")
    axes[-1].set_xlabel("Window start position")
    fig.suptitle("Real dataset: pi along chromosome")
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    fig.savefig(FIG_SUPP / "supplementary_real_pi_along_chromosome.png", dpi=300)
    plt.close(fig)

    # Call R script for additional figures
    r_script = ROOT / "scripts" / "make_figures.R"
    if r_script.exists():
        print("Generating R-based figures...")
        result = subprocess.run(["Rscript", str(r_script)], check=False)
        if result.returncode != 0:
            print("Warning: R script execution failed, but continuing with Python figures", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
