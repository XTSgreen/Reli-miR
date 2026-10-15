# E4 multi-cohort transfer figure.
# Panel A: Spearman rho (internal Rel vs external per-miRNA Pearson) per cohort,
#          point estimate + B=2000 patient-bootstrap CI; all positive.
# Panel B: external median Pearson for all shared miRNAs vs Rel>=0.7 top-33 subset.
# Run from the repository root:  Rscript new/code/gen_fig_E4_multicohort.R
source("new/code/figure_style.R")

s <- read.csv("new/results/E4_multicohort/multicohort_E4_summary.csv",
              stringsAsFactors = FALSE, check.names = FALSE)
s <- s[order(s$rho_Pearson_point_sample), ]

export_plot("F_E4_multicohort_transfer", function() {
  par(mfrow = c(1, 2), mar = c(4.5, 9.5, 1.6, .8), mgp = c(2.9, .7, 0),
      cex = .85, las = 1)
  y <- seq_len(nrow(s))
  plot(NA, xlim = c(0, .68), ylim = c(.4, nrow(s) + .6),
       yaxt = "n", xlab = "Spearman rho: internal Rel vs external Pearson",
       ylab = "", bty = "l")
  axis(2, at = y, labels = paste0(s$cohort, " (n=", s$n_patients, ")"), las = 1)
  abline(v = 0, lty = 2, col = "#888888")
  arrows(s$rho_Pearson_low, y, s$rho_Pearson_high, y,
         angle = 90, length = .04, code = 3, col = "#555555")
  points(s$rho_Pearson_point_sample, y, pch = 16, cex = 1.1, col = pal["Rel"])
  panel_letter("A")

  plot(NA, xlim = c(.4, nrow(s) + .6), ylim = c(0, .68),
       xaxt = "n", xlab = "", ylab = "External median Pearson", bty = "l")
  axis(1, at = y, labels = s$cohort, las = 1)
  for (i in y) {
    segments(i - .18, s$overall_median_Pearson[i], i + .18, s$overall_median_Pearson[i],
             col = "#777777", lwd = 2.4)
    points(i, s$top33_median_Pearson[i], pch = 17, cex = 1.15, col = pal["Rel"])
    segments(i, s$overall_median_Pearson[i], i, s$top33_median_Pearson[i], lty = 3, col = "#BBBBBB")
  }
  legend("topright", c("Rel>=0.7 (33 items)", "All shared miRNAs"),
         col = c(pal["Rel"], "#777777"), pch = c(17, NA), lty = c(NA, 1),
         lwd = c(NA, 2.4), bty = "n", cex = .9)
  panel_letter("B")
}, width = 9.5, height = 5.2)
