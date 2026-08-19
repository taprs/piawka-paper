#!/usr/bin/env Rscript
# Figure 3: piawka-vs-pixy Hudson's Fst percentiles on the real North/South
# A. lyrata dataset, for two gene sets (GWAS candidate genes, GO:2000028
# flowering-time genes). Reads the rank tables written by
# scripts/analyze_northsouth_fst_ranks.py. One point per gene per tool (per
# gene-interval, for the rare genes with >1 genomic interval in genes.bed),
# sorted left to right from highest to lowest piawka percentile.
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
fig_main <- file.path(root, "figures", "main")
dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw())

# Per-tool percentile dotplot: piawka and pixy percentiles for each gene,
# sorted by piawka percentile, colored by tool. Returns the plot plus the gene
# count, so the two panels can be width-weighted in the combined figure.
percentile_plot <- function(table_name, label) {
  d <- read_tsv(file.path(tables_dir, table_name), show_col_types = FALSE) %>%
    mutate(
      gene_label = ifelse(!is.na(common_name) & common_name != "", common_name, gene_lyrata),
      gene_label = make.unique(gene_label),
      gene_label = factor(gene_label, levels = gene_label[order(-piawka_percentile)])
    )
  d_long <- d %>%
    select(gene_label, piawka_percentile, pixy_percentile) %>%
    pivot_longer(c(piawka_percentile, pixy_percentile), names_to = "tool", values_to = "percentile") %>%
    mutate(tool = sub("_percentile$", "", tool))
  list(
    plot = ggplot(d_long, aes(x = gene_label, y = percentile, color = tool)) +
      geom_point(size = 2, alpha = 0.7, shape = 16) +
      scale_color_manual(values = scales::hue_pal()(2)) +
      labs(
        x = NULL,
        y = expression(paste(F[ST], " percentile, north vs south individuals")),
        color = "Tool",
        title = label
      ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, face = "italic")),
    n = nrow(d)
  )
}

gwas_pct <- percentile_plot("northsouth_fst_ranks_gwas_snp_genes.tsv", "GWAS genes")
go_pct <- percentile_plot(
  "northsouth_fst_ranks_go2000028_genes.tsv",
  "GO:2000028 genes (regulation of photoperiodism, flowering)"
)

# Figure 3: both panels side by side (GWAS left, GO right), shared legend.
# Panel widths are proportional to gene counts so the x-axis spacing matches.
f3 <- (gwas_pct$plot + go_pct$plot) +
  plot_layout(guides = "collect", axes = "collect", widths = c(gwas_pct$n, go_pct$n)) &
  theme(legend.position = "bottom")
ggsave(file.path(fig_main, "f3.png"), f3, width = 8, height = 4, dpi = 600, limitsize = FALSE)

cat("Figure 3 generated successfully\n")
