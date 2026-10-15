# Identity-and-expression figure (revised manuscript figure 10).
# Panels: (a) reliability vs input expression mean; (b) reliability vs input variance;
#         (c) host-gene co-expression by matched stratum (Rel >= 0.7 vs background);
#         (d) stratified label-permutation null for the observed difference.
# The permutation reproduces biological_annotation.R exactly (seed 20261014, B = 10000).
#
# Reads only published result tables. Set RELI_TABLES_ROOT to the tables directory
# (default: results_tables).
style <- if (file.exists("new/code/figure_style.R")) "new/code/figure_style.R" else "figures/figure_style.R"
source(style)
tables_root <- Sys.getenv("RELI_TABLES_ROOT", "results_tables")
if (!dir.exists(tables_root)) stop("tables root not found: ", tables_root, "; set RELI_TABLES_ROOT")
tb <- function(...) file.path(tables_root, ...)

z <- read.csv(tb("E5", "E5_miRNA_identity_expression_score.csv"), check.names = FALSE)
assoc <- read.csv(tb("E5", "E5_host_expression_associations_all_mRNA.csv"), check.names = FALSE)
published <- read.csv(tb("E5", "E5_host_expression_stratified_label_permutation.csv"), check.names = FALSE)
stopifnot(nrow(z) == 846)

matched <- z[is.finite(z$Rel) & is.finite(z$native_expression_mean_366) &
             is.finite(z$native_expression_variance_366), , drop = FALSE]
matched$score_group <- ifelse(matched$Rel >= .7, "high", "background")
matched$mean_bin <- pmin(10L, ceiling(rank(matched$native_expression_mean_366,
                                           ties.method = "average") / nrow(matched) * 10))
matched$variance_bin <- pmin(10L, ceiling(rank(log1p(matched$native_expression_variance_366),
                                              ties.method = "average") / nrow(matched) * 10))
matched$stratum <- paste(matched$mean_bin, matched$variance_bin, sep = "_")
assoc_abs <- aggregate(abs(host_expression_correlation) ~ miRNA, assoc, mean)
names(assoc_abs)[2] <- "mean_abs_host_rho"
matched <- merge(matched, assoc_abs, by = "miRNA", all.x = TRUE)
high <- matched[matched$score_group == "high" & is.finite(matched$mean_abs_host_rho), , drop = FALSE]
backg <- matched[matched$score_group == "background" & is.finite(matched$mean_abs_host_rho), , drop = FALSE]
strata <- sort(intersect(unique(high$stratum), unique(backg$stratum)))
common <- matched[matched$stratum %in% strata & is.finite(matched$mean_abs_host_rho), , drop = FALSE]
nh <- sum(common$score_group == "high")
weights <- vapply(strata, function(st) sum(common$stratum == st & common$score_group == "high") / nh, numeric(1))
per_stratum <- t(vapply(strata, function(st) {
  s <- common[common$stratum == st, , drop = FALSE]
  c(high = mean(s$mean_abs_host_rho[s$score_group == "high"]),
    background = mean(s$mean_abs_host_rho[s$score_group == "background"]),
    n_high = sum(s$score_group == "high"), n_bg = sum(s$score_group == "background"))
}, numeric(4)))
observed <- sum(weights * (per_stratum[, "high"] - per_stratum[, "background"]))

B <- 10000L
set.seed(20261014L)
null <- rep(NA_real_, B)
for (b in seq_len(B)) {
  null[b] <- sum(weights * vapply(strata, function(st) {
    s <- common[common$stratum == st, , drop = FALSE]
    k <- sum(s$score_group == "high")
    idx <- sample.int(nrow(s), size = k, replace = FALSE)
    mean(s$mean_abs_host_rho[idx]) - mean(s$mean_abs_host_rho[-idx])
  }, numeric(1)))
}
p_rep <- (1 + sum(abs(null) >= abs(observed))) / (B + 1)
pub <- published[published$analysis_scope == "all_18204_mRNA", ]
cat(sprintf("reproduced difference = %.7f (published %.7f)\n", observed,
            pub$observed_stratum_weighted_mean_abs_rho_difference))
cat(sprintf("reproduced p = %.6f (published %.6f); strata = %d\n", p_rep,
            pub$label_permutation_p_two_sided, length(strata)))

export_plot("F_E5_identity_features", function() {
  par(mfrow = c(2, 2), mar = c(4.3, 4.6, 1.8, .7), mgp = c(2.9, .7, 0), cex = .85, las = 1)

  ## (a)/(b) reliability against input expression scale
  for (v in c("native_expression_mean_366", "native_expression_variance_366")) {
    h <- z$Rel >= .7
    plot(log1p(z[[v]]), z$Rel, pch = ifelse(h, 17, 16), cex = ifelse(h, .8, .45),
         col = ifelse(h, pal[["Pearson"]], adjustcolor(pal[["Rel"]], .3)),
         xlab = if (v == "native_expression_mean_366") "log(1 + input expression mean)" else
                 "log(1 + input expression variance)",
         ylab = "Reliability (Rel)", ylim = c(0, 1), bty = "l")
    legend("topleft", c("Rel >= 0.7 (n=33)", "Other miRNAs (n=813)"), pch = c(17, 16),
           col = c(pal[["Pearson"]], pal[["Rel"]]), bty = "n", cex = .78)
    panel_letter(if (v == "native_expression_mean_366") "A" else "B")
  }

  ## (c) host-gene co-expression by matched stratum
  yy <- seq_along(strata)
  rng <- range(per_stratum[, c("high", "background")])
  plot(NA, xlim = rng + c(-.02, .02), ylim = c(.4, length(strata) + .6), yaxt = "n",
       xlab = "Mean |Spearman rho| of host-gene pair", ylab = "", bty = "l")
  axis(2, at = yy, labels = sprintf("%s (h=%d/b=%d)", strata, per_stratum[, "n_high"],
       per_stratum[, "n_bg"]), las = 1, cex.axis = .68)
  for (i in yy) {
    segments(per_stratum[i, "background"], i, per_stratum[i, "high"], i, col = "#CCCCCC", lwd = 1.6)
    points(per_stratum[i, "background"], i, pch = 1, cex = .95, col = "#777777")
    points(per_stratum[i, "high"], i, pch = 16, cex = .95, col = pal[["Rel"]])
  }
  legend("topleft", c("Rel >= 0.7", "Background"), col = c(pal[["Rel"]], "#777777"),
         pch = c(16, 1), bty = "n", cex = .78)
  panel_letter("C")

  ## (d) stratified label-permutation null
  hist(null, breaks = 60, col = "#DDDDDD", border = "white", freq = FALSE, main = "",
       xlab = "Stratum-weighted difference in mean |rho|", ylab = "Density", bty = "l")
  abline(v = observed, col = pal[["Rel"]], lwd = 2.2)
  legend("topleft", c(sprintf("Observed %.4f", observed),
                      sprintf("p = %.3f (%s permutations)", p_rep, format(B, big.mark = ","))),
         col = c(pal[["Rel"]], NA), lwd = c(2.2, NA), bty = "n", cex = .78)
  panel_letter("D")
}, width = 9.6, height = 7.6)

cat("figure written:", file.path(figure_root, "F_E5_identity_features.pdf"), "\n")
