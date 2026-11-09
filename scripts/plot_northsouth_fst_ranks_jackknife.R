#!/usr/bin/env Rscript
# Figure 3 variant (percentile panels A, B only) with 95% CI bars from the leave-one-South-out (7-fold) jackknife
# written by scripts/analyze_a3_jackknife.py. Genes whose piawka and pixy CIs do not overlap are highlighted.
# Writes figures/main/f3.png (compact version: all GWAS genes, GO genes with piawka percentile >= 75) and figures/supp/S3_flowering_time_fst.png (all GO:2000028 genes, no GWAS panel).
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(script_dir) == 0) {
  script_dir <- getwd()
}
root <- dirname(script_dir)
tables_dir <- file.path(root, "results", "tables")
fig_supp <- file.path(root, "figures", "supp")
fig_main <- file.path(root, "figures", "main")
dir.create(fig_supp, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw())

percentile_plot <- function(table_name, label, min_piawka_pct = -Inf, sc = 1, sites = FALSE) {
  # sc scales point size, bar width and axis-label size with the canvas
  # sites = TRUE overlays stacked bars of per-gene site counts on a secondary y-axis: pixy sites used (pixy_snps) with
  # the extra sites only piawka uses (piawka_snps - pixy_snps) stacked on top
  d <- read_tsv(file.path(tables_dir, table_name), show_col_types = FALSE) %>%
    mutate(
      gene_label = ifelse(!is.na(common_name) & common_name != "", common_name, gene_lyrata),
      gene_label = make.unique(gene_label),
      gene_label = factor(gene_label, levels = gene_label[order(-piawka_percentile)])
    )
  # the table is restricted to piawka percentile >= min_piawka_pct; the highlighted genes are chosen from the full set
  # highlighted genes: those whose piawka and pixy jackknife CIs do not overlap
  top <- d %>%
    filter(piawka_ci_lo > pixy_ci_hi | pixy_ci_lo > piawka_ci_hi)
  d <- d %>%
    filter(piawka_percentile >= min_piawka_pct) %>%
    mutate(gene_label = droplevels(gene_label))
  top <- filter(top, gene_label %in% d$gene_label)
  d_long <- d %>%
    select(gene_label, piawka_percentile, pixy_percentile, ends_with("_ci_lo"), ends_with("_ci_hi")) %>%
    pivot_longer(-gene_label, names_to = c("tool", ".value"), names_pattern = "(piawka|pixy)_(percentile|ci_lo|ci_hi)")
  k <- max(d$piawka_snps) / 100  # sites per percentile unit on the shared axis
  bars <- d %>%
    transmute(gene_label, `pixy sites used` = pixy_snps / k, `piawka extra sites` = (piawka_snps - pixy_snps) / k) %>%
    pivot_longer(-gene_label, names_to = "part", values_to = "y") %>%
    mutate(part = factor(part, levels = c("pixy sites used", "piawka extra sites")))
  site_layers <- if (sites) list(
    geom_col(data = bars, aes(x = gene_label, y = y, fill = part), inherit.aes = FALSE,
             width = 0.7, alpha = 0.45, position = position_stack(reverse = TRUE)),
    scale_fill_manual(values = c("pixy sites used" = "grey35", "piawka extra sites" = "grey80"), name = "Sites"),
    scale_y_continuous(sec.axis = sec_axis(~ . * k, name = "Sites used"))
  ) else list()
  # highlighted genes get a dashed vertical line beneath their points and a bold x-axis label
  # both tools' bars sit on the gene's own gridline, equally thick and semi-transparent so overlaps show as a blend
  list(
    plot = ggplot(d_long, aes(x = gene_label, y = percentile, color = tool)) +
      site_layers +
      geom_segment(data = top, aes(x = gene_label, xend = gene_label, y = -Inf, yend = Inf),
                   inherit.aes = FALSE, color = "grey45", linewidth = 0.5, linetype = "dashed") +
      geom_linerange(data = filter(d_long, tool == "piawka"), aes(ymin = ci_lo, ymax = ci_hi),
                     linewidth = 0.8 * sc, alpha = 0.65) +
      geom_linerange(data = filter(d_long, tool == "pixy"), aes(ymin = ci_lo, ymax = ci_hi),
                     linewidth = 0.8 * sc, alpha = 0.65) +
      geom_point(size = 1.6 * sc, alpha = 0.7, shape = 16) +
      scale_x_discrete(limits = levels(d$gene_label)) +
      expand_limits(y = c(0, 100)) +
      scale_color_manual(values = scales::hue_pal()(2)) +
      labs(
        x = NULL,
        y = expression(paste(F[ST], " percentile, north vs south individuals")),
        color = "Tool",
        title = label
      ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6 * sc, face = ifelse(levels(d$gene_label) %in% top$gene_label, "bold.italic", "italic"))),
    n = nrow(d)
  )
}

make_fig <- function(min_piawka_pct, out_path, width = 8, height = 4, base_size = 11) {
  sc <- base_size / 11
  theme_set(theme_bw(base_size = base_size))
  gwas_pct <- percentile_plot("northsouth_jackknife_ci_gwas_snp_genes.tsv", "Latitude GWAS genes", sc = sc)  # all GWAS genes, never filtered
  go_pct <- percentile_plot(
    "northsouth_jackknife_ci_go2000028_genes.tsv",
    # narrower panel in the filtered version: wrap the title so it is not clipped
    ifelse(is.finite(min_piawka_pct), "Flowering time genes\n(>75th percentile)",
           "GO:2000028 genes (regulation of photoperiodism, flowering)"),
    min_piawka_pct = min_piawka_pct, sc = sc
  )
  if (is.finite(min_piawka_pct)) {
    # short canvas: break the y-axis title over two lines so it is not clipped
    ylab2 <- expression(atop(F[ST] ~ "percentile,", "north vs south individuals"))
    gwas_pct$plot <- gwas_pct$plot + labs(y = ylab2)
    go_pct$plot <- go_pct$plot + labs(y = ylab2)
  }
  # same layout and size as Figure 3: both panels side by side, shared legend, widths proportional to gene counts
  top_row <- (gwas_pct$plot + go_pct$plot) +
    plot_layout(guides = "collect", axes = "collect", widths = c(gwas_pct$n, go_pct$n)) &
    theme(legend.position = "bottom", legend.key.height = unit(0.8 * sc, "lines"))
  if (is.finite(min_piawka_pct)) top_row <- top_row & theme(panel.grid.minor.y = element_blank(), legend.margin = margin(0, 0, 0, 0))
  ggsave(out_path, wrap_elements(top_row), width = width, height = height, dpi = 600, limitsize = FALSE)
  cat(basename(out_path), ": ", gwas_pct$n, " GWAS genes, ", go_pct$n, " GO genes\n", sep = "")
}

# supplementary: all GO:2000028 genes only (the GWAS genes are shown in full in Figure 3), with the site-count overlay
go_all <- percentile_plot("northsouth_jackknife_ci_go2000028_genes.tsv", "GO:2000028 genes (regulation of photoperiodism, flowering)", sites = TRUE)
ggsave(file.path(fig_supp, "S3_flowering_time_fst.png"), go_all$plot + theme(legend.position = "bottom"),
       width = 8, height = 4, dpi = 600, limitsize = FALSE)
cat("S3_flowering_time_fst.png: ", go_all$n, " GO genes\n", sep = "")
# the compact version (all GWAS genes; GO genes with piawka >= 75th percentile) is the paper's Figure 3
make_fig(75, file.path(fig_main, "f3.png"), width = 5.5, height = 3)
