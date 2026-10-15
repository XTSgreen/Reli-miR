# E1 part 2 (fixed): Jaccard overlap, and per-split downstream (E3) comparison.
out <- "new/results/E1_evaluation"
seeds <- 20261014:20261033
models <- c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")
pm <- "median_patient_correlation"
cm <- "correlation_matrix_mean_absolute_difference"

cat("=== F. Jaccard overlap of selected sets between splits ===\n")
jr <- read.csv(file.path(out, "E1_Rel_threshold_pairwise_Jaccard.csv"), check.names = FALSE)
for (m in models) {
  for (th in c(0.5, 0.6, 0.7)) {
    v <- jr$Jaccard[jr$model == m & jr$threshold == th]
    cat(sprintf("Rel>=%.1f %-15s Jaccard: median %.3f, range [%.3f, %.3f]\n",
      th, m, median(v), min(v), max(v)))
  }
}
je <- read.csv(file.path(out, "E1_pairwise_Jaccard.csv"), check.names = FALSE)
je <- je[je$method == "Rel", ]
cat("-- equal-size top-3.9%/8.7%/20.9% --\n")
for (rate in unique(je$coverage_rate)) {
  for (m in models) {
    v <- je$Jaccard[je$coverage_rate == rate & je$model == m]
    cat(sprintf("  coverage %.3f %-15s Jaccard median %.3f [%.3f, %.3f]\n",
      rate, m, median(v), min(v), max(v)))
  }
}

cat("\n=== G/H. per-split downstream (E3) vs random and vs simple scores ===\n")
e3 <- do.call(rbind, lapply(seeds, function(s) {
  p <- file.path(out, "by_seed", as.character(s), "E3_winner_equal_size_summary.csv")
  if (!file.exists(p)) return(NULL)
  d <- read.csv(p, check.names = FALSE); d$seed <- s
  d[, c("seed", "subset", "n_miRNA", pm, cm,
        "median_variance_ratio", "DE_effect_slope_pred_vs_real", "Cox_logHR_correlation")]
}))
e3$is_random <- grepl("^Random_", e3$subset)
cat("rows:", nrow(e3), "| seeds:", length(unique(e3$seed)), "| sizes:", paste(sort(unique(e3$n_miRNA)), collapse = ","), "\n\n")

for (k in c(33, 74, 177)) {
  d <- e3[e3$n_miRNA == k, ]
  cat(sprintf("--- size %d (coverage %.3f) ---\n", k, unique(d$requested_coverage)))
  # per-seed random mean
  rand <- aggregate(cbind(pm_r = d[[pm]][d$is_random], cm_r = d[[cm]][d$is_random]),
                    list(seed = d$seed[d$is_random]), mean)
  named <- d[!d$is_random, ]
  for (s in c("Rel", "Pearson", "calRMSE", "training_OOF_RMSE")) {
    v <- named[named$subset == s, ]
    cmp <- merge(v[, c("seed", pm, cm)], rand, by = "seed")
    cat(sprintf("  %-18s patient-corr %.4f (mean over seeds; random %.4f) [above random in %d/%d] | CMAD %.4f (random %.4f) [below in %d/%d]\n",
      s, mean(cmp[[pm]]), mean(cmp$pm_r),
      sum(cmp[[pm]] > cmp$pm_r, na.rm = TRUE), nrow(cmp),
      mean(cmp[[cm]]), mean(cmp$cm_r),
      sum(cmp[[cm]] < cmp$cm_r, na.rm = TRUE), nrow(cmp)))
  }
  cat(sprintf("  %-18s patient-corr %.4f | CMAD %.4f (100-200 random lists per split)\n",
    "Random(mean)", mean(rand$pm_r), mean(rand$cm_r)))
  # Rel vs Pearson head-to-head per seed
  wide <- reshape(named[, c("seed", "subset", pm, cm)], idvar = "seed",
                  timevar = "subset", direction = "wide")
  if (all(c(paste0(pm, ".Rel"), paste0(pm, ".Pearson")) %in% names(wide))) {
    dd <- wide[[paste0(pm, ".Rel")]] - wide[[paste0(pm, ".Pearson")]]
    dc <- wide[[paste0(cm, ".Rel")]] - wide[[paste0(cm, ".Pearson")]]
    cat(sprintf("  Rel-Pearson per-seed: patient-corr diff mean %+.4f (Rel higher in %d/%d); CMAD diff mean %+.4f (Rel lower in %d/%d)\n",
      mean(dd, na.rm = TRUE), sum(dd > 0, na.rm = TRUE), sum(is.finite(dd)),
      mean(dc, na.rm = TRUE), sum(dc < 0, na.rm = TRUE), sum(is.finite(dc))))
  }
  # amplitude shrinkage in Rel subset
  rv <- named[named$subset == "Rel", ]
  cat(sprintf("  Rel subset: median variance ratio %.4f; DE slope %.4f; Cox logHR corr %.4f\n",
    mean(rv$median_variance_ratio, na.rm = TRUE),
    mean(rv$DE_effect_slope_pred_vs_real, na.rm = TRUE),
    mean(rv$Cox_logHR_correlation, na.rm = TRUE)))
  cat("\n")
}
