#!/usr/bin/env Rscript
# Statistical tests on real-data biallelic piawka vs pixy per-window pi
# (results/tables/data_gain_pi_vs_site_gain.tsv):
#
#   1. Correlation between sites gained by piawka over pixy (absolute and
#      percent) and the corresponding pi gain (piawka - pixy), per window.
#   2. Equality of the coefficients of variation of per-window pi between
#      piawka and pixy (Feltz & Miller asymptotic test; Krishnamoorthy & Lee
#      modified signed-likelihood ratio test), via the cvequality package.
#
# Output: results/tables/real_pi_site_gain_stats_tests.tsv

suppressMessages({
  library(cvequality)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(script_dir) == 0) {
  script_dir <- getwd()
}
root <- dirname(script_dir)
tables_dir <- file.path(root, "results", "tables")

df <- read.delim(file.path(tables_dir, "data_gain_pi_vs_site_gain.tsv"))
df$pixy_pi <- as.numeric(df$pixy_pi)
df$piawka_pi <- as.numeric(df$piawka_pi)
df$pi_gain_piawka_vs_pixy <- as.numeric(df$pi_gain_piawka_vs_pixy)
df$sites_gained_piawka_vs_pixy <- as.numeric(df$sites_gained_piawka_vs_pixy)
df$sites_gained_pct_vs_pixy <- as.numeric(df$sites_gained_pct_vs_pixy)

results <- list()

# 1. Sites gained (absolute, %) vs pi gain -----------------------------------
sub_abs <- df[complete.cases(df[, c("pi_gain_piawka_vs_pixy", "sites_gained_piawka_vs_pixy")]), ]
p_abs <- cor.test(sub_abs$sites_gained_piawka_vs_pixy, sub_abs$pi_gain_piawka_vs_pixy, method = "pearson")
s_abs <- cor.test(sub_abs$sites_gained_piawka_vs_pixy, sub_abs$pi_gain_piawka_vs_pixy, method = "spearman")

sub_pct <- df[complete.cases(df[, c("pi_gain_piawka_vs_pixy", "sites_gained_pct_vs_pixy")]), ]
p_pct <- cor.test(sub_pct$sites_gained_pct_vs_pixy, sub_pct$pi_gain_piawka_vs_pixy, method = "pearson")
s_pct <- cor.test(sub_pct$sites_gained_pct_vs_pixy, sub_pct$pi_gain_piawka_vs_pixy, method = "spearman")

results[[length(results) + 1]] <- data.frame(
  test = "pearson", comparison = "sites_gained_abs_vs_pi_gain",
  n = nrow(sub_abs), statistic = unname(p_abs$estimate), p_value = p_abs$p.value
)
results[[length(results) + 1]] <- data.frame(
  test = "spearman", comparison = "sites_gained_abs_vs_pi_gain",
  n = nrow(sub_abs), statistic = unname(s_abs$estimate), p_value = s_abs$p.value
)
results[[length(results) + 1]] <- data.frame(
  test = "pearson", comparison = "sites_gained_pct_vs_pi_gain",
  n = nrow(sub_pct), statistic = unname(p_pct$estimate), p_value = p_pct$p.value
)
results[[length(results) + 1]] <- data.frame(
  test = "spearman", comparison = "sites_gained_pct_vs_pi_gain",
  n = nrow(sub_pct), statistic = unname(s_pct$estimate), p_value = s_pct$p.value
)

# 2. Equality of coefficient of variation of per-window pi -------------------
sub_cv <- df[complete.cases(df[, c("pixy_pi", "piawka_pi")]) & df$pixy_pi > 0 & df$piawka_pi > 0, ]
vals <- c(sub_cv$pixy_pi, sub_cv$piawka_pi)
grp <- factor(c(rep("pixy", nrow(sub_cv)), rep("piawka", nrow(sub_cv))))

asymp <- asymptotic_test(x = vals, y = grp)
mslr <- mslr_test(nr = 1e4, x = vals, y = grp)

# The CV values themselves, as quoted in the paper text (percent).
cv_pct <- function(x) 100 * sd(x) / mean(x)
results[[length(results) + 1]] <- data.frame(
  test = "cv_percent", comparison = "pi_cv_piawka",
  n = nrow(sub_cv), statistic = cv_pct(sub_cv$piawka_pi), p_value = NA_real_
)
results[[length(results) + 1]] <- data.frame(
  test = "cv_percent", comparison = "pi_cv_pixy",
  n = nrow(sub_cv), statistic = cv_pct(sub_cv$pixy_pi), p_value = NA_real_
)

results[[length(results) + 1]] <- data.frame(
  test = "cv_equality_asymptotic", comparison = "pi_cv_piawka_vs_pixy",
  n = nrow(sub_cv), statistic = asymp$D_AD, p_value = asymp$p_value
)
results[[length(results) + 1]] <- data.frame(
  test = "cv_equality_mslr", comparison = "pi_cv_piawka_vs_pixy",
  n = nrow(sub_cv), statistic = mslr$MSLRT, p_value = mslr$p_value
)

out <- do.call(rbind, results)
out_path <- file.path(tables_dir, "real_pi_site_gain_stats_tests.tsv")
write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote", out_path, "\n")
print(out)
