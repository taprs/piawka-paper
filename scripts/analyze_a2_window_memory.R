#!/usr/bin/env Rscript
# A2: does peak memory (max RSS) scale linearly with window size? Simulated window series, 1 process.
# Per tool: linear fit (slope, 95% CI, R2), Pearson/Spearman correlation, log-log slope, lack-of-fit F-test
# (linear vs window size as a factor). Output: results/tables/a2_window_memory_fit.tsv
library(dplyr)
library(readr)

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(script_dir) == 0) script_dir <- getwd()
root <- dirname(script_dir)
d <- read_tsv(file.path(root, "results", "tables", "a2_scaling.tsv"), show_col_types = FALSE) %>%
  filter(series == "sim_windows")

fit <- function(x) {
  m <- lm(max_rss_mb ~ window_mbp, x)
  ci <- confint(m)["window_mbp", ]
  ll <- lm(log(max_rss_mb) ~ log(window_mbp), x)
  llci <- confint(ll)["log(window_mbp)", ]
  pe <- cor.test(x$window_mbp, x$max_rss_mb)
  sp <- suppressWarnings(cor.test(x$window_mbp, x$max_rss_mb, method = "spearman"))
  tibble(n = nrow(x), rss_at_min_mb = mean(x$max_rss_mb[x$window_mbp == min(x$window_mbp)]),
         rss_at_max_mb = mean(x$max_rss_mb[x$window_mbp == max(x$window_mbp)]),
         intercept_mb = coef(m)[[1]], slope_mb_per_mbp = coef(m)[[2]], slope_ci_lo = ci[[1]], slope_ci_hi = ci[[2]],
         r2 = summary(m)$r.squared, slope_p = coef(summary(m))[2, 4],
         pearson_r = pe$estimate, pearson_p = pe$p.value, spearman_rho = sp$estimate, spearman_p = sp$p.value,
         loglog_slope = coef(ll)[[2]], loglog_ci_lo = llci[[1]], loglog_ci_hi = llci[[2]],
         lack_of_fit_p = anova(m, lm(max_rss_mb ~ factor(window_mbp), x))$`Pr(>F)`[2])
}
res <- d %>% group_by(tool) %>% group_modify(~ fit(.x)) %>% ungroup()
write_tsv(res, file.path(root, "results", "tables", "a2_window_memory_fit.tsv"), na = "")
print(as.data.frame(res), digits = 4)
