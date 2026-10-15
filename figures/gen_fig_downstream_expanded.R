# Downstream preservation and equal-size controls (revised manuscript figure 8).
# Panels: (a) fidelity metrics across reliability thresholds;
#         (b) equal-size control on within-patient correlation;
#         (c) equal-size control on correlation-matrix mean absolute difference;
#         (d) amplitude diagnostics across subset sizes (mean over the 20 repeated splits).
#
# Reads only published result tables. Set RELI_TABLES_ROOT to the tables directory
# (default: results_tables).
style <- if (file.exists("new/code/figure_style.R")) "new/code/figure_style.R" else "figures/figure_style.R"
source(style)
tables_root <- Sys.getenv("RELI_TABLES_ROOT", "results_tables")
if (!dir.exists(tables_root)) stop("tables root not found: ", tables_root, "; set RELI_TABLES_ROOT")
tb <- function(...) file.path(tables_root, ...)

ss <- read.csv(tb("E3", "E3_subset_summary_named.csv"), check.names = FALSE)
rd <- read.csv(tb("E3", "E3_random_subset_metrics_slim.csv"), check.names = FALSE)
ps <- read.csv(tb("E3", "E3_downstream_by_split.csv"), check.names = FALSE)
named <- function(k, s) ss[ss$n_miRNA == k & ss$subset == s, ][1, ]

sets <- list(All = named(846, "All"), R05 = named(177, "Rel"),
             R06 = named(74, "Rel"), R07 = named(33, "Rel"))
metrics <- c(median_patient_correlation = "Median within-patient correlation",
             correlation_matrix_correlation = "Correlation-matrix preservation",
             DE_effect_correlation = "Differential-effect correlation",
             Cox_logHR_correlation = "Cox logHR correlation")
mcol <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
sizes <- c(33, 74, 177)
subsets <- c("Rel", "Pearson", "RMSE")

## values used in panel (a) and panel (d), computed outside the plotting closure for reporting
y <- vapply(names(metrics), function(m)
  vapply(sets, function(d) as.numeric(d[[m]])[1], numeric(1)), numeric(length(sets)))
rownames(y) <- names(sets); colnames(y) <- unname(metrics)
stopifnot(all(is.finite(y)))
rows_ps <- ps[ps$subset == "Rel", ]
vr <- sapply(sizes, function(k) mean(rows_ps$median_variance_ratio[rows_ps$n_miRNA == k], na.rm = TRUE))
sl <- sapply(sizes, function(k) mean(rows_ps$DE_effect_slope_pred_vs_real[rows_ps$n_miRNA == k], na.rm = TRUE))

export_plot("F_downstream_preservation_expanded", function() {
  par(mfrow = c(2, 2), mar = c(4.3, 4.6, 1.8, .7), mgp = c(2.9, .7, 0), cex = .85, las = 1)

  ## (a) fidelity metrics across thresholds
  matplot(1:4, y, type = "b", lty = 1, pch = c(16, 17, 15, 18), col = mcol,
          xaxt = "n", xlab = "Reliability threshold (retained 846 / 177 / 74 / 33)",
          ylab = "Preservation", ylim = c(.3, 1), bty = "l")
  axis(1, at = 1:4, labels = c("All", "0.5", "0.6", "0.7"), cex.axis = .9)
  legend("bottomright", colnames(y), col = mcol, pch = c(16, 17, 15, 18), lty = 1,
         bty = "n", cex = .78)
  panel_letter("A")

  ## (b)/(c) equal-size controls
  for (mi in list(c("median_patient_correlation", "Median within-patient correlation", "B"),
                  c("correlation_matrix_mean_absolute_difference",
                    "Correlation-matrix mean absolute difference", "C"))) {
    met <- mi[1]; lab <- mi[2]; letter <- mi[3]
    rand <- lapply(sizes, function(k) rd[[met]][rd$requested_size == k])
    ylim <- if (met == "median_patient_correlation") c(.1, .92) else c(.1, .32)
    boxplot(rand, names = sizes, ylim = ylim, outline = FALSE, col = "#EEEEEE",
            border = "#888888", xlab = "Number of selected miRNAs", ylab = lab, bty = "l")
    off <- c(Rel = -.08, Pearson = 0, RMSE = .08)
    for (s in subsets) {
      v <- vapply(sizes, function(k) as.numeric(named(k, s)[[met]])[1], numeric(1))
      points(seq_along(sizes) + off[s], v, col = pal[s],
             pch = c(Rel = 16, Pearson = 17, RMSE = 15)[s], cex = 1.15)
    }
    legend("bottomright", c("Rel", "Pearson", "RMSE", "Random (1,000 lists)"),
           col = c(pal[["Rel"]], pal[["Pearson"]], pal[["RMSE"]], "#888888"),
           pch = c(16, 17, 15, 0), bty = "n", cex = .8)
    panel_letter(letter)
  }

  ## (d) amplitude diagnostics, averaged over the 20 repeated splits
  plot(sizes, vr, type = "b", pch = 16, col = "#0072B2", ylim = c(0, 1),
       xlab = "Number of selected miRNAs", ylab = "Value", bty = "l")
  lines(sizes, sl, type = "b", pch = 17, col = "#D55E00")
  abline(h = 1, lty = 3, col = "#888888")
  legend("topright", c("Median variance ratio", "Differential-effect slope"),
         col = c("#0072B2", "#D55E00"), pch = c(16, 17), lty = 1, bty = "n", cex = .82)
  panel_letter("D")
}, width = 9.6, height = 7.6)

cat("figure written:", file.path(figure_root, "F_downstream_preservation_expanded.pdf"), "\n")
cat("panel a:\n"); print(round(y, 4))
cat("panel d variance ratio (20-split mean):", round(vr, 4), "\n")
cat("panel d effect slope (20-split mean):", round(sl, 4), "\n")
