source("new/code/figure_style.R")
z <- read_result("E2/E2_Rel_threshold_coverage_FNR_point.csv")
b <- read_result("E2/E2_FNR_bootstrap_summary.csv")
b <- b[abs(b$Pearson_label_threshold - .5) < 1e-9, ]
stopifnot(all(b$B_valid_FNR == 2000), nrow(z) == nrow(b))
export_plot("F_E2_quality_coverage", function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.3, 1.4, .8), mgp = c(2.7, .7, 0), cex = .9, las = 1)
  plot(z$threshold, z$precision, type = "o", pch = 16, col = pal["Rel"], ylim = c(0, 1), xlab = "Rel threshold", ylab = "Fraction", bty = "l")
  lines(z$threshold, z$recall, type = "o", pch = 17, col = pal["Pearson"])
  lines(z$threshold, z$coverage, type = "o", pch = 1, lty = 2, col = pal["Random"])
  legend("bottomleft", c("Precision", "Recall", "Coverage"), col = c(pal["Rel"], pal["Pearson"], pal["Random"]), pch = c(16, 17, 1), lty = c(1, 1, 2), bty = "n")
  panel_letter("A")
  plot(z$threshold, z$FNR, type = "n", ylim = c(0, 1), xlab = "Rel threshold", ylab = "Predictive false-negative fraction", bty = "l")
  polygon(c(b$threshold, rev(b$threshold)), c(b$FNR_lo, rev(b$FNR_hi)), col = adjustcolor(pal["Rel"], .18), border = NA)
  lines(z$threshold, z$FNR, type = "o", pch = 16, col = pal["Rel"])
  for (th in c(.5, .6, .7)) {
    a <- z[which.min(abs(z$threshold - th)), ]
    text(a$threshold + .045, a$FNR - .055, paste0("m = ", a$n_selected), cex = .78)
  }
  legend("bottomright", c("Original-sample estimate", "Joint bootstrap 95% interval"), lty = c(1, NA), pch = c(16, 15), col = c(pal["Rel"], adjustcolor(pal["Rel"], .35)), bty = "n", cex = .95)
  panel_letter("B")
}, width = 8.5, height = 3.8)
