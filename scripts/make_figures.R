#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(stringr)

# Determine root and output paths using the script's location
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(script_dir) == 0) {
  script_dir <- getwd()
}
root <- dirname(script_dir)
tables_dir <- file.path(root, "results", "tables")
fig_main <- file.path(root, "figures", "main")
dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw())

# Figure 1: resources, 16 rescaled neutral one-population replicates (scripts/run_resource_sim16.sh -> analyze_a2_scaling.py)
scal <- read_tsv(file.path(tables_dir, "a2_scaling.tsv"), show_col_types = FALSE)

# Figure 1A: peak memory vs window size (10 Mbp analysed as 0.625-10 Mbp windows, 1 process); thin lines = replicates, thick = median
fw <- scal %>% filter(series == "sim_windows")
f1a <- ggplot(fw, aes(window_mbp, max_rss_mb, color = tool)) +
  geom_line(aes(group = interaction(replicate, tool)), alpha = 0.2) +
  geom_point(alpha = 0.4, shape = 1) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2.5) +
  scale_x_log10(breaks = c(0.625, 1.25, 2.5, 5, 10), labels = c("0.625", "1.25", "2.5", "5", "10")) +
  scale_y_log10() +
  labs(x = "Window size, Mbp", y = "Max RAM usage, MB", color = "Tool")

# Figure 1B: runtime vs number of parallel processes (eight 1.25 Mbp windows)
ft <- scal %>% filter(series == "sim_threads")
f1b <- ggplot(ft, aes(factor(processes), elapsed_sec, color = tool)) +
  geom_line(aes(group = interaction(replicate, tool)), alpha = 0.2) +
  geom_point(alpha = 0.4, shape = 1) +
  stat_summary(fun = median, geom = "line", aes(group = tool), linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2.5) +
  labs(x = "# parallel processes", y = "Seconds elapsed", color = "Tool")

# Figure 1 pileup
f1_fin <- f1a + f1b + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")
ggsave(file.path(fig_main, "f1.png"), f1_fin, width = 8, height = 3.5, dpi = 600)

# Figure 2: data gain + pi gain given REF/ALT-agnostic filtering
f2m <- read_tsv(file.path(tables_dir, "multiallelic_pi_real_pairs.tsv"))

f2ab <- read_tsv(file.path(tables_dir, "data_gain_pi_vs_site_gain.tsv")) %>%
  filter(!grepl("MIX|mix", group))

f2 <- left_join(f2ab, f2m) %>%
  pivot_longer(matches("_pi$"), names_to = "tool", values_to = "pi") %>%
  pivot_longer(contains("_sites"), names_to = "tool2", values_to = "sites") %>%
  filter(
    substr(tool, 1, 5) == substr(tool2, 1, 5),
    tool != "bial_pi",
    !grepl("MIX|mix", group)
  ) %>%
  mutate(tool = str_extract(tool, '^[^_]+'))

# Figure 2A: pixy vs piawka sites
f2a <- ggplot(f2ab, aes(pixy_sites, piawka_sites, color = reorder(group, -piawka_pi))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "Nucleotides used / 10Kbp\n(pixy)",
    y = "Nucleotides used / 10Kbp\n(piawka)",
    color = "Population"
  ) +
  theme(legend.position = 'none')

# Figure 2B: pixy vs piawka pi
f2b <- ggplot(f2ab, aes(pixy_pi, piawka_pi, color = reorder(group, -piawka_pi))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10(labels = function(x) sprintf("%.2f", x * 100)) +
  scale_y_log10(labels = function(x) sprintf("%.2f", x * 100)) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "% nucleotide diversity / 10Kbp\n(pixy)",
    y = "% nucleotide diversity / 10Kbp\n(piawka)",
    color = "Population"
  ) +
  theme(legend.position = 'none')

# Figure 2C: sites vs pi
f2c <- ggplot(f2, aes(sites, pi, color = reorder(group, -pi), shape = tool)) +
  geom_line(aes(group = interaction(chromosome, window_pos_1, group)), alpha=.03) +
  geom_point(alpha = 0.2) +
  geom_point(
    aes(group = group), size = 4, stroke = 1.2,
    data = f2 %>% group_by(group, tool) %>% summarize(sites = median(sites), pi = median(pi), .groups = "drop")
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(
    limits = c("piawka", "mult", "pixy", "pixymult"),
    labels = c("piawka, biallelic", "piawka, multiallelic", "pixy, biallelic", "pixy, multiallelic"),
    values = c(2, 6, 0, 5)
  ) +
  scale_y_log10(labels = function(x) sprintf("%.2f", x * 100)) +
  guides(shape = guide_legend(order = 2), color = guide_legend(order = 1)) +
  labs(
    x = "Nucleotides used / 10Kbp\n(polymorphic + invariant)",
    y = "% nucleotide diversity / 10Kbp",
    color = "Population",
    shape = "Tool, mode"
  )

# Figure 2 pileup
f2_fin <- (f2a + f2b) / f2c + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect", heights = c(2,3))
ggsave(file.path(fig_main, "f2.png"), f2_fin, width = 8, height = 7, dpi = 600)

cat("R figures generated successfully\n")
