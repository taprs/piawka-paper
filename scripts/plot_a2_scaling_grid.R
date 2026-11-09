#!/usr/bin/env Rscript
# Supplementary Figure S1 (A2): resource usage vs samples / populations / window size, real and simulated data.
# Rows (top to bottom): real memory, real time, simulated memory, simulated time.
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(script_dir) == 0) script_dir <- getwd()
root <- dirname(script_dir)
d <- read_tsv(file.path(root, "results", "tables", "a2_scaling.tsv"), show_col_types = FALSE) %>%
  mutate(procs = factor(processes, c(1, 8), c("1 process", "8 processes")))
# same colours as Figure 1 (ggplot default hue palette: piawka, pixy)
cols <- setNames(scales::hue_pal()(2), c("piawka", "pixy"))
theme_set(theme_bw())

med <- function(x, ...) x %>% group_by(tool, procs, ...) %>%
  summarise(across(c(elapsed_sec, max_rss_mb), median), .groups = "drop")

lab <- list(
  rss = "Peak memory (MB, log scale)",
  time = "Wall time (s)"
)

panel <- function(dat, x, y, xlab, ylab, logx = FALSE, logy = FALSE, breaks = waiver(), shape = NULL, smooth = FALSE) {
  if (is.null(dat) || nrow(dat) == 0) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = "not measured", colour = "grey40") +
             labs(x = xlab, y = ylab) + theme(axis.text = element_blank(), axis.ticks = element_blank()))
  }
  g <- ggplot(dat, aes(.data[[x]], .data[[y]], colour = tool, linetype = procs)) +
    geom_line() + geom_point() + scale_colour_manual(values = cols) +
    scale_linetype_manual(values = c("1 process" = "solid", "8 processes" = "dashed"),
                          limits = c("1 process", "8 processes")) +
    labs(x = xlab, y = ylab, colour = NULL, linetype = NULL)
  # linear regression of RSS on window size (fit on untransformed values), one line per tool x process count.
  # Axes are log-scaled via coord_trans (not scale_*_log10) so geom_smooth fits in linear space; the lines look curved on log axes.
  if (smooth) {
    g <- g + geom_smooth(aes(group = interaction(tool, procs)), method = "lm", formula = y ~ x, se = TRUE,
                         colour = "grey30", fill = "grey50", alpha = 0.25, linewidth = 0.6, linetype = "dotted") +
      scale_x_continuous(breaks = breaks, labels = function(b) format(b, trim = TRUE, scientific = FALSE, drop0trailing = TRUE)) +
      coord_trans(x = "log10", y = "log10")
    return(g)
  }
  if (!is.null(shape)) g <- g + aes(shape = .data[[shape]], group = interaction(tool, procs, .data[[shape]])) + labs(shape = NULL)
  g <- if (logx) g + scale_x_log10(breaks = breaks, labels = function(b) format(b, trim = TRUE, scientific = FALSE, drop0trailing = TRUE)) else g + scale_x_continuous(breaks = breaks)
  if (logy) g + scale_y_log10() else g + expand_limits(y = 0)
}

real_s <- med(filter(d, series == "real_samples"), n_samples)
real_p <- med(filter(d, series == "real_groups"), n_groups)
sim_s <- med(filter(d, series == "sim_small"), series, n_samples) %>%
  mutate(series = factor(series, "sim_small", "100-sample replicates"))
sim_w <- med(filter(d, series == "sim_windows"), window_mbp)

# rows: real memory, real time, simulated memory, simulated time; columns: samples, populations (real) / window size (simulated)
real_row <- function(y, ylab, logy) list(
  panel(real_s, "n_samples", y, "Samples (log scale; 12 populations)", ylab, logx = TRUE, logy = logy),
  panel(real_p, "n_groups", y, "Populations (log scale; 910 samples)", ylab, logx = TRUE, logy = logy, breaks = c(1, 2, 4, 8, 12)))
sim_row <- function(y, ylab, logy) list(
  panel(sim_s, "n_samples", y, "Samples (log scale)", ylab, logx = TRUE, logy = logy, shape = "series"),
  panel(sim_w, "window_mbp", y, "Window size (Mbp, log scale; 100 samples)", ylab, logx = TRUE, logy = logy,
        breaks = c(0.625, 1.25, 2.5, 5, 10)))
ps <- c(real_row("max_rss_mb", lab$rss, TRUE), real_row("elapsed_sec", lab$time, FALSE),
        sim_row("max_rss_mb", lab$rss, TRUE), sim_row("elapsed_sec", lab$time, FALSE))
# collect one legend; the simulated-samples panel (5th) also carries the 100-sample-replicates key
ps[-5] <- lapply(ps[-5], function(g) g + theme(legend.position = "none"))
p <- wrap_plots(ps, ncol = 2) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
dir.create(file.path(root, "figures", "supp"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(root, "figures", "supp", "S1_resource_usage.png"), p, width = 8.5, height = 13, dpi = 200)
