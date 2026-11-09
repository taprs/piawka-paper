#!/usr/bin/env Rscript
# Supplementary figure (A4): Figure 2C data as ratios, per population (the populations of Figure 2).
#  A: pi on multiallelic / biallelic sites within each tool (sina plots; piawka/pixy colours as in Figures 1 and 3)
#  B: piawka / pixy pi on biallelic sites (default mode)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(script_dir) == 0) script_dir <- getwd()
root <- dirname(script_dir)
theme_set(theme_bw())

f2ab <- read_tsv(file.path(root, "results", "tables", "data_gain_pi_vs_site_gain.tsv"), show_col_types = FALSE) %>%
  filter(!grepl("MIX|mix", group))
pop_order <- f2ab %>% group_by(group) %>% summarise(m = mean(piawka_pi)) %>% arrange(-m) %>% pull(group)  # as Fig 2A/B
pop_col <- setNames(RColorBrewer::brewer.pal(8, "Dark2")[seq_along(pop_order)], pop_order)  # Fig 2 colours
pop_lab <- sprintf("<span style='color:%s'>\u25cf</span> %s", pop_col, pop_order)

d <- read_tsv(file.path(root, "results", "tables", "multiallelic_pi_real_pairs.tsv"), show_col_types = FALSE) %>%
  filter(group %in% pop_order, bial_pi > 0, pixy_pi > 0, mult_pi > 0, pixymult_pi > 0) %>%
  mutate(group = factor(group, pop_order), xi = as.integer(group))

a <- d %>% transmute(group, xi, piawka = mult_pi / bial_pi, pixy = pixymult_pi / pixy_pi) %>%
  pivot_longer(c(piawka, pixy), names_to = "tool", values_to = "ratio") %>%
  mutate(tool = factor(tool, c("piawka", "pixy")), x = xi + ifelse(tool == "piawka", -0.23, 0.23))
b <- d %>% transmute(group, xi, x = xi, ratio = bial_pi / pixy_pi)

write_tsv(bind_rows(a %>% mutate(panel = "multiallelic/biallelic") %>% rename(label = tool),
                    b %>% mutate(panel = "piawka/pixy, biallelic", label = "piawka/pixy")) %>%
            group_by(panel, group, label) %>%
            summarise(median_ratio = median(ratio), q25 = quantile(ratio, .25), q75 = quantile(ratio, .75), .groups = "drop"),
          file.path(root, "results", "tables", "multiallelic_pi_ratio_by_population.tsv"))

xscale <- scale_x_continuous(breaks = seq_along(pop_order), labels = pop_lab, expand = expansion(add = 0.6))
med <- stat_summary(fun = median, geom = "crossbar", width = 0.4, linewidth = 0.3, colour = "black", fatten = 1)
ytick <- scale_y_continuous(labels = function(x) sprintf("%.1f", x))
yline <- geom_hline(yintercept = 1, linewidth = 0.6, colour = "black")

pa <- ggplot(a, aes(x, ratio, colour = tool, group = interaction(group, tool))) +
  yline + ggforce::geom_sina(maxwidth = 0.46, scale = FALSE, bw = 0.03, size = 0.4, alpha = 0.4) + med +
  scale_colour_manual(values = scales::hue_pal()(2), guide = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  ytick + xscale + labs(x = NULL, y = "π multiallelic / π biallelic", colour = "Tool") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

pb <- ggplot(b, aes(x, ratio, group = group)) +
  yline +
  ggforce::geom_sina(maxwidth = 0.8, scale = FALSE, bw = 0.15, size = 0.4, alpha = 0.4, colour = "grey40") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, linewidth = 0.3, colour = "black", fatten = 1) +
  ytick + xscale + labs(x = NULL, y = "π piawka / π pixy\n(biallelic sites)") +
  theme(axis.text.x = ggtext::element_markdown(angle = 40, hjust = 1))

p <- pa / pb + plot_annotation(tag_levels = "A") + plot_layout(heights = c(1, 1))
dir.create(file.path(root, "figures", "supp"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(root, "figures", "supp", "S2_pi_ratios.png"), p, width = 7, height = 7, dpi = 300)
