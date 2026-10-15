source("new/code/figure_style.R")
z <- read_result("E5/E5_miRNA_identity_expression_score.csv")
ex <- read_result("E5/E5_high_mid_low_example_scatter_data.csv")
stopifnot(nrow(z) == 846, nrow(ex) == 432)
export_plot("F_E5_expression_associations", function() {
  par(mfrow = c(1, 2), mar = c(4.3, 4.3, 1.4, .8), mgp = c(2.7, .7, 0), cex = .9, las = 1)
  for (v in c("native_expression_mean_366", "native_expression_variance_366")) {
    high <- z$Rel >= .7
    plot(log1p(z[[v]]), z$Rel, pch = ifelse(high, 17, 16), cex = ifelse(high, .8, .5), col = ifelse(high, pal["Pearson"], adjustcolor(pal["Rel"], .35)), ylab = "Calibration Rel", xlab = if (v == "native_expression_mean_366") "log(1 + input expression mean)" else "log(1 + input expression variance)", bty = "l", ylim = c(0, 1))
    legend("topleft", c("Rel >= 0.7 (n=33)", "Other miRNAs (n=813)"), pch = c(17, 16), col = c(pal["Pearson"], pal["Rel"]), bty = "n", cex = .9)
    panel_letter(if (v == "native_expression_mean_366") "A" else "B")
  }
}, width = 8.5, height = 3.8)
export_plot("F_E5_prediction_examples", function() {
  par(mfrow = c(1, 3), mar = c(4.4, 4.1, 2.1, .6), mgp = c(2.5, .7, 0), cex = 1, las = 1)
  for (i in seq_along(unique(ex$miRNA))) {
    a <- ex[ex$miRNA == unique(ex$miRNA)[i], ]
    lim <- range(a$observed_native, a$predicted_native, finite = TRUE)
    plot(a$observed_native, a$predicted_native, xlim = lim, ylim = lim, pch = ifelse(a$split == "test", 16, 1), col = ifelse(a$split == "test", pal["Rel"], "#777777"), cex = .6, xlab = "Observed input expression", ylab = "Predicted input expression", bty = "l")
    abline(a = 0, b = 1, lty = 2, col = "#999999")
    mtext(paste0(LETTERS[i], "  ", unique(a$miRNA), "; Rel = ", round(unique(a$Rel), 3)), side = 3, adj = 0, line = .4, cex = .95)
    legend("topleft", c("Calibration (n=72)", "Test (n=72)"), pch = c(1, 16), col = c("#777777", pal["Rel"]), bty = "n", cex = .95)
  }
}, width = 10, height = 3.6)
