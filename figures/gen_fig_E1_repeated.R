# E1 repeated-split stability figure.
# Panel A: paired model differences in median test Pearson across seeds.
# Panel B: retained miRNA counts at Rel>=0.5/0.6/0.7 across seeds, per model.
# Panel C: selection frequency of the original 33 high-reliability miRNAs
#          into each split's top-3.9% Rel list, per model.
# Run from the repository root:  Rscript new/code/gen_fig_E1_repeated.R
source("new/code/figure_style.R")

perf <- read_result("E1_evaluation/E1_model_performance_by_seed.csv")
orig33 <- read_result("E1_evaluation/E1_original_33_selection_frequency.csv")

seed_dirs <- list.dirs("new/results/E1_evaluation/by_seed", recursive = FALSE)
thr <- do.call(rbind, lapply(seed_dirs, function(d) {
  p <- file.path(d, "Rel_threshold_counts.csv")
  if (!file.exists(p)) return(NULL)
  x <- read.csv(p, check.names = FALSE)
  x[x$selection_kind == "Rel_threshold", ]
}))
n_seeds <- length(unique(perf$seed))
message("E1 figure: seeds with evaluation output = ", n_seeds)

models <- c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")
model_pal <- c(ElasticNet = "#0072B2", PCA_ElasticNet = "#D55E00", PCA_Ridge = "#009E73")
model_pch <- c(ElasticNet = 16, PCA_ElasticNet = 17, PCA_Ridge = 15)

export_plot("F_E1_repeated_stability", function() {
  par(mfrow = c(1, 3), mar = c(4.5, 4.4, 1.6, .6), mgp = c(2.9, .7, 0),
      cex = .9, las = 1)

  ## Panel A: paired model differences (median test Pearson)
  contrasts <- list(
    c("PCA_ElasticNet", "PCA_Ridge"),
    c("ElasticNet", "PCA_ElasticNet"),
    c("ElasticNet", "PCA_Ridge"))
  diffs <- lapply(contrasts, function(p) {
    a <- perf[perf$seed %in% perf$seed & perf$model == p[1], c("seed", "median_pearson")]
    b <- perf[perf$model == p[2], c("seed", "median_pearson")]
    m <- merge(a, b, by = "seed", suffixes = c("_a", "_b"))
    m$median_pearson_a - m$median_pearson_b
  })
  ylim <- range(unlist(diffs), na.rm = TRUE)
  ylim <- ylim + c(-.06, .06) * diff(ylim)
  plot(NA, xlim = c(.5, 3.5), ylim = ylim, xaxt = "n", xlab = "", ylab = "Delta median test Pearson",
       bty = "l")
  axis(1, at = 1:3, labels = c("PCAEN\n- PCAR", "EN\n- PCAEN", "EN\n- PCAR"), tick = TRUE, las = 1)
  abline(h = 0, lty = 2, col = "#888888")
  for (i in seq_along(diffs)) {
    points(rep(i, length(diffs[[i]])) + jitter(rep(0, length(diffs[[i]])), amount = .05),
           diffs[[i]], pch = 16, cex = .8, col = "#333333")
    points(i, mean(diffs[[i]], na.rm = TRUE), pch = 4, cex = 1.6, col = "#000000", lwd = 2)
  }
  panel_letter("A")

  ## Panel B: retained counts per Rel threshold, per model
  thr_levels <- c(0.5, 0.6, 0.7)
  plot(NA, xlim = c(.4, 3.6), ylim = range(thr$n_selected, na.rm = TRUE) * c(0, 1.15),
       xaxt = "n", xlab = "Rel threshold", ylab = "Retained miRNAs", bty = "l")
  axis(1, at = 1:3, labels = c("0.5", "0.6", "0.7"))
  for (mi in seq_along(models)) {
    m <- models[mi]
    off <- c(-.22, 0, .22)[mi]
    for (ti in seq_along(thr_levels)) {
      v <- thr$n_selected[thr$model == m & thr$threshold == thr_levels[ti]]
      points(rep(ti + off, length(v)) + jitter(rep(0, length(v)), amount = .04), v,
             pch = model_pch[mi], col = adjustcolor(model_pal[mi], .55), cex = .75)
      segments(ti + off - .1, median(v, na.rm = TRUE), ti + off + .1, median(v, na.rm = TRUE),
               col = model_pal[mi], lwd = 2.4)
    }
  }
  legend("topright", models, col = model_pal, pch = model_pch, bty = "n", cex = .85)
  panel_letter("B")

  ## Panel C: original 33 miRNAs' frequency in each split's top-3.9% Rel list
  freq_col <- "Rel_top_0039_frequency"
  plot(NA, xlim = c(-.02, 1.02), ylim = c(.4, 3.6),
       xlab = "Selection frequency into top-3.9% Rel list", yaxt = "n", ylab = "", bty = "l")
  axis(2, at = seq_along(models), labels = c("EN", "PCAEN", "PCAR"), las = 1)
  abline(v = c(0, 1), lty = 3, col = "#BBBBBB")
  med <- numeric(length(models))
  for (mi in seq_along(models)) {
    m <- models[mi]
    v <- orig33[[freq_col]][orig33$model == m]
    v <- v[is.finite(v)]
    med[mi] <- median(v, na.rm = TRUE)
    points(v, rep(mi, length(v)) + jitter(rep(0, length(v)), amount = .13),
           pch = model_pch[mi], col = adjustcolor(model_pal[mi], .55), cex = .8)
    points(med[mi], mi, pch = 4, cex = 1.5, col = "#000000", lwd = 2)
  }
  panel_letter("C")
}, width = 11.5, height = 4)
