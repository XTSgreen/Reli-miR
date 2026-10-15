# Calibration-and-threshold figure (revised manuscript figure 4).
# Panels: (a) per-miRNA reliability vs test RMSE with threshold lines;
#         (b) score ablation with joint-resampling intervals;
#         (c) precision / recall / coverage across thresholds;
#         (d) predicted false-negative rate with resampling intervals.
#
# Reads only published result tables. Set RELI_TABLES_ROOT to the tables directory
# (default: results_tables).
#   Rscript verification/...   ->   Rscript figures/gen_fig_E2_calibration_thresholds.R
style <- if (file.exists("new/code/figure_style.R")) "new/code/figure_style.R" else "figures/figure_style.R"
source(style)
tables_root <- Sys.getenv("RELI_TABLES_ROOT", "results_tables")
if (!dir.exists(tables_root)) stop("tables root not found: ", tables_root, "; set RELI_TABLES_ROOT")
tb <- function(...) file.path(tables_root, ...)

pm <- read.csv(tb("E2", "E2_per_miRNA_reliability_test_error.csv"), check.names = FALSE)
stopifnot(nrow(pm) == 846)
rho_cal <- cor(pm$Rel, pm$test_RMSE_trainingSD, method = "spearman")
cat(sprintf("panel a: Spearman(Rel, test RMSE) = %.7f (manuscript reports -0.5978)\n", rho_cal))

sb <- read.csv(tb("E2", "E2_complete_score_joint_bootstrap_summary.csv"), check.names = FALSE)
cv <- read.csv(tb("E2", "E2_Rel_threshold_coverage_FNR_point.csv"), check.names = FALSE)
fn <- read.csv(tb("E2", "E2_FNR_bootstrap_summary.csv"), check.names = FALSE)
fn <- fn[abs(fn$Pearson_label_threshold - .5) < 1e-9, ]

export_plot("F_E2_calibration_thresholds", function() {
  par(mfrow = c(2, 2), mar = c(4.3, 4.6, 1.8, .7), mgp = c(2.9, .7, 0), cex = .85, las = 1)

  ## (a) reliability vs test error
  plot(pm$test_RMSE_trainingSD, pm$Rel, pch = 16, cex = .42,
       col = adjustcolor(pal[["Rel"]], .3),
       xlab = "Test RMSE (training-SD scale)", ylab = "Reliability (Rel)", bty = "l")
  abline(h = c(.5, .6, .7), lty = 3, col = "#BBBBBB")
  text(max(pm$test_RMSE_trainingSD), c(.5, .6, .7), c("0.5", "0.6", "0.7"),
       adj = c(1, -.4), cex = .72, col = "#888888")
  legend("topright", sprintf("Spearman = %.4f", rho_cal), bty = "n", cex = .82)
  panel_letter("A")

  ## (b) score ablation with joint-resampling intervals
  par(mar = c(4.3, 8.8, 1.8, .5))
  d <- sb[order(sb$point_rho_vs_test_RMSE), ]
  yy <- seq_len(nrow(d))
  xr <- range(c(d$rho_CI_q025, d$rho_CI_q975), na.rm = TRUE)
  plot(NA, xlim = c(xr[1] - .09, xr[2] + .01), ylim = c(.4, nrow(d) + .6), yaxt = "n",
       xlab = "Spearman vs test RMSE (resampling interval)", ylab = "", bty = "l")
  axis(2, at = yy, labels = d$score, las = 1, cex.axis = .62)
  abline(v = 0, lty = 3, col = "#BBBBBB")
  arrows(d$rho_CI_q025, yy, d$rho_CI_q975, yy, angle = 90, length = .03, code = 3, col = "#888888")
  points(d$point_rho_vs_test_RMSE, yy, pch = 16, cex = 1.0, col = pal[["Rel"]])
  points(d$bootstrap_median_rho, yy, pch = 0, cex = .9, col = "#333333")
  legend("topleft", c("Observed", "Resampling median", "95% interval"),
         pch = c(16, 0, NA), lty = c(NA, NA, 1), col = c(pal[["Rel"]], "#333333", "#888888"),
         bty = "n", cex = .68)
  panel_letter("B")

  ## (c) precision / recall / coverage across thresholds
  par(mar = c(4.3, 4.6, 1.8, .7))
  plot(NA, xlim = c(0, .92), ylim = c(0, 1), xlab = "Reliability threshold",
       ylab = "Value", bty = "l")
  cols <- c(precision = "#D55E00", recall = "#009E73", coverage = "#0072B2")
  for (v in names(cols)) {
    keep <- is.finite(cv[[v]])
    lines(cv$threshold[keep], cv[[v]][keep], type = "b", pch = 16, cex = .5, col = cols[v])
  }
  legend("right", c("Precision", "Recall", "Coverage"), col = cols, pch = 16, lty = 1,
         bty = "n", cex = .78)
  panel_letter("C")

  ## (d) predicted false-negative rate with intervals
  keep <- is.finite(fn$FNR_lo) & is.finite(fn$FNR_hi)
  plot(NA, xlim = c(0, .92), ylim = c(0, 1), xlab = "Reliability threshold",
       ylab = "Predicted false-negative rate", bty = "l")
  polygon(c(fn$threshold[keep], rev(fn$threshold[keep])),
          c(fn$FNR_lo[keep], rev(fn$FNR_hi[keep])), col = "#E8F0F8", border = NA)
  lines(fn$threshold[keep], fn$FNR[keep], type = "b", pch = 16, cex = .5, col = pal[["Rel"]])
  legend("bottomright", c("Predicted FNR", "95% resampling interval"),
         col = c(pal[["Rel"]], "#E8F0F8"), pch = c(16, NA), lty = c(1, NA),
         fill = c(NA, "#E8F0F8"), border = NA, bty = "n", cex = .78)
  panel_letter("D")
}, width = 9.6, height = 7.6)

cat("figure written:", file.path(figure_root, "F_E2_calibration_thresholds.pdf"), "\n")
