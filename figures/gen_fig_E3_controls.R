source("new/code/figure_style.R")
z <- read_result("E3/E3_subset_summary.csv")
sizes <- c(33, 74, 177)
export_plot("F_E3_equal_size_controls", function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.4, 1.4, .8), mgp = c(2.9, .7, 0), cex = .9, las = 1)
  for (metric in c("median_patient_correlation", "correlation_matrix_mean_absolute_difference")) {
    rand <- lapply(sizes, function(m) z[z$requested_size == m & grepl("^Random", z$subset), metric])
    stopifnot(all(lengths(rand) == 1000))
    ylim <- if (metric == "median_patient_correlation") c(.1, .92) else c(.07, .43)
    boxplot(rand, names = sizes, ylim = ylim, outline = FALSE, col = "#EEEEEE", border = "#888888", xlab = "Number of selected miRNAs", ylab = if (metric == "median_patient_correlation") "Median within-patient correlation" else "Mean absolute correlation-matrix difference", bty = "l")
    for (s in c("Rel", "Pearson", "RMSE")) {
      a <- z[match(paste(s, sizes), paste(z$subset, z$requested_size)), metric]
      points(seq_along(sizes) + c(Rel = -.08, Pearson = 0, RMSE = .08)[s], a, col = pal[s], pch = c(Rel = 16, Pearson = 17, RMSE = 15)[s], cex = 1.1)
    }
    legend(if (metric == "median_patient_correlation") "topright" else "bottomleft", c("Rel", "Pearson", "RMSE", "Random (1,000 lists)"), col = pal, pch = c(16, 17, 15, 0), bty = "n", cex = .95)
    panel_letter(if (metric == "median_patient_correlation") "A" else "B")
  }
}, width = 9, height = 4)
