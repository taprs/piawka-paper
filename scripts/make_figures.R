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
fig_supp <- file.path(root, "figures", "supplementary")
dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_supp, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw())

# Figure 1: resources
f1 <- read_tsv(file.path(tables_dir, 'resource_usage_timing.tsv')) %>%
  mutate(
    scenario = ifelse(substr(sample_id, 15, 20) == "bneckL", "bottleneck", "neutral"),
    scenario = factor(scenario, c("neutral", "bottleneck"))
  )

fm <- read_tsv(file.path(tables_dir, 'resource_usage_1thread_synwin_vs_threads8.tsv')) %>%
  pivot_longer(contains("max_rss_kb")) %>%
  mutate(
    window_size = ifelse(name == "max_rss_kb_syn_win", "10Mbp", "1.25Mbp"),
    scenario = ifelse(substr(sample_id, 15, 20) == "bneckL", "bottleneck", "neutral"),
    scenario = factor(scenario, c("neutral", "bottleneck"))
  )

# Figure 1A: memory use @ two window sizes
f1a <- ggplot(fm, aes(factor(window_size), value, color = tool)) +
  geom_line(aes(group = interaction(sample_id, tool, scenario), linetype = scenario), alpha = 0.2) +
  geom_point(aes(shape = scenario)) +
  scale_y_log10(labels = function(x) x / 1e6) +
  scale_shape_manual(values = c(1, 4)) +
  labs(x = "Window size, Mbp", y = "Max RAM usage, GB", color = "Tool")

# Figure 1B: runtime
f1b <- ggplot(f1, aes(factor(threads), elapsed_sec, color = tool)) +
  geom_line(aes(group = interaction(sample_id, tool, scenario), linetype = scenario), alpha = 0.2) +
  geom_point(aes(shape = scenario)) +
  scale_shape_manual(values = c(1, 4)) +
  labs(x = "# parallel processes", y = "Seconds elapsed", color = "Tool")

# Figure 1 pileup
f1_fin <- f1a + f1b + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")
ggsave(file.path(fig_main, "f1.png"), f1_fin, width = 8, height = 4, dpi = 300)

# Supplementary figure 1: memory use at different #threads
ggplot(f1, aes(factor(threads), max_rss_kb, color = tool)) +
  geom_point(position = 'identity') +
  geom_line(aes(group = interaction(sample_id, tool))) +
  facet_grid(cols = vars(scenario)) +
  scale_y_log10(labels = function(x) x / 1e6) +
  labs(x = "# parallel processes", y = "Max RAM usage, GB", color = "Tool")
ggsave(file.path(fig_supp, "supplementary_resource_threads.png"), width = 10, height = 4, dpi = 300)

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
  geom_point(alpha = 0.2) +
  geom_point(
    aes(group = group), size = 4, stroke = 1.2,
    data = f2 %>% group_by(group, tool) %>% summarize(sites = median(sites), pi = median(pi), .groups = "drop")
  ) +
  scale_shape_manual(
    limits = c("pixy", "piawka", "mult"),
    labels = c("pixy, biallelic", "piawka, biallelic", "piawka, multiallelic"),
    values = c(0, 1, 2)
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_log10(labels = function(x) sprintf("%.2f", x * 100)) +
  labs(
    x = "Nucleotides used / 10Kbp\n(polymorphic + invariant)",
    y = "% nucleotide diversity / 10Kbp",
    color = "Population",
    shape = "Tool, mode"
  )

# Figure 2 pileup
f2_fin <- (f2a + f2b) / f2c + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")
ggsave(file.path(fig_main, "f2.png"), f2_fin, width = 6, height = 4, dpi = 300)

# Figure 3: π-thetaW vs π-thetaLow
f3 <- read_tsv(file.path(tables_dir, "accuracy_theta_proxy_per_run.tsv")) %>%
  select(!pixy_pi_minus_theta_w) %>%
  pivot_longer(matches("minus_.*w$")) %>%
  mutate(scenario = ifelse(source_class == "bneck", "bottleneck", "neutral"))

# Figure 3A: thetaW
f3a <- ggplot(filter(f3, name == "piawka_pi_minus_theta_w"), aes(factor(missing_rate * 100), value, color = reorder(scenario, -value))) +
  geom_point() +
  geom_line(aes(group = interaction(sample_id, name))) +
  labs(x = "% missing data", y = expression(pi - theta[Watterson]), color = "Simulation scenario") +
  scale_color_brewer(palette = "Set1", direction = -1)

# Figure 3B: thetaLow
f3b <- ggplot(filter(f3, name == "piawka_pi_minus_theta_low"), aes(factor(missing_rate * 100), value, color = reorder(scenario, -value))) +
  geom_point() +
  geom_line(aes(group = interaction(sample_id, name))) +
  labs(x = "% missing data", y = expression(pi - theta[low]), color = "Simulation scenario") +
  scale_color_brewer(palette = "Set1", direction = -1)

# Figure 3 pileup
f3_fin <- f3a / f3b + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect", axes = "collect")
ggsave(file.path(fig_main, "f3.png"), f3_fin, width = 5, height = 5, dpi = 300)

cat("R figures generated successfully\n")
