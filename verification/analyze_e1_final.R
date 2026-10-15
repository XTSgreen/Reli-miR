# Comprehensive E1 analysis for the revision report.
# 20 pre-registered patient-level splits (seeds 20261014-20261033); splits overlap,
# so no ordinary t-tests across splits: report distributions and per-seed directions.
out <- "new/results/E1_evaluation"
perf <- read.csv(file.path(out, "E1_model_performance_by_seed.csv"), check.names = FALSE)
qual <- read.csv(file.path(out, "E1_selection_quality_by_seed.csv"), check.names = FALSE)
seeds <- sort(unique(perf$seed))
n_seed <- length(seeds)
models <- c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")

cat("=== A. model ranking stability (20 splits) ===\n")
# Winner by each criterion
for (crit in c("median_pearson", "mean_r2", "normalized_rmse")) {
  w <- sapply(seeds, function(s) {
    d <- perf[perf$seed == s, ]
    if (crit == "normalized_rmse") d$model[which.min(d[[crit]])] else d$model[which.max(d[[crit]])]
  })
  tb <- table(factor(w, levels = models))
  cat(sprintf("%-16s winner counts: %s\n", crit,
    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", ")))
}
cat("training-internal selected model counts: ",
    paste(sprintf("%s=%d", names(table(factor(perf$training_selected_model, levels = models))),
                  as.integer(table(factor(perf$training_selected_model, levels = models)))), collapse = ", "), "\n")

# Paired per-seed differences (test metrics)
cat("\n--- per-seed paired differences (mean of 20, and #seeds favoring first) ---\n")
contrasts <- list(c("ElasticNet", "PCA_ElasticNet"), c("ElasticNet", "PCA_Ridge"),
                  c("PCA_ElasticNet", "PCA_Ridge"))
for (crit in c("median_pearson", "mean_r2", "normalized_rmse")) {
  for (cc in contrasts) {
    dl <- sapply(seeds, function(s) {
      a <- perf$seed == s & perf$model == cc[1]; b <- perf$seed == s & perf$model == cc[2]
      v <- perf[[crit]][a] - perf[[crit]][b]
      if (crit == "normalized_rmse") -v else v
    })
    cat(sprintf("%-16s %s-%s: mean diff %+.4f; first better in %d/%d splits\n", crit,
      cc[1], cc[2], mean(dl), sum(dl > 0), n_seed))
  }
}
cat("\n--- per-model distribution summary (median [q05, q95]) ---\n")
for (m in models) {
  d <- perf[perf$model == m, ]
  cat(sprintf("%-15s median_pearson %.4f [%.4f, %.4f]; mean_r2 %.4f [%.4f, %.4f]; norm_rmse %.4f [%.4f, %.4f]\n",
    m, mean(d$median_pearson), quantile(d$median_pearson, .05), quantile(d$median_pearson, .95),
    mean(d$mean_r2), quantile(d$mean_r2, .05), quantile(d$mean_r2, .95),
    mean(d$normalized_rmse), quantile(d$normalized_rmse, .05), quantile(d$normalized_rmse, .95)))
}

cat("\n=== B. reliability score vs test quality (per-seed Spearman) ===\n")
assoc <- do.call(rbind, lapply(seeds, function(s) {
  p <- file.path(out, "by_seed", as.character(s), "score_test_associations.csv")
  if (!file.exists(p)) return(NULL)
  read.csv(p, check.names = FALSE)
}))
for (m in models) {
  d <- assoc[assoc$model == m, ]
  cat(sprintf("\n-- %s --\n", m))
  methods <- unique(d$method)
  for (mm in methods) {
    v <- d$spearman_quality_score_vs_test_RMSE_trainingSD[d$method == mm]
    v <- v[is.finite(v)]
    if (!length(v)) next
    cat(sprintf("  %-28s rho vs RMSE: mean %+.4f, median %+.4f, range [%+.4f, %+.4f], negative in %d/%d\n",
      mm, mean(v), median(v), min(v), max(v), sum(v < 0), length(v)))
  }
}

cat("\n=== C. retained counts at Rel thresholds ===\n")
thr <- do.call(rbind, lapply(seeds, function(s) {
  p <- file.path(out, "by_seed", as.character(s), "Rel_threshold_counts.csv")
  if (!file.exists(p)) return(NULL)
  read.csv(p, check.names = FALSE)
}))
thr <- thr[thr$selection_kind == "Rel_threshold", ]
for (m in models) {
  for (t in c(0.5, 0.6, 0.7)) {
    v <- thr$n_selected[thr$model == m & thr$threshold == t]
    cov <- thr$actual_coverage[thr$model == m & thr$threshold == t]
    cat(sprintf("%-15s Rel>=%.1f: n median %5.0f range [%4.0f, %4.0f]; coverage median %.4f\n",
      m, t, median(v), min(v), max(v), median(cov)))
  }
}

cat("\n=== D. original-33 membership stability (aggregated over 20 splits) ===\n")
o33 <- read.csv(file.path(out, "E1_original_33_selection_frequency.csv"), check.names = FALSE)
for (m in models) {
  d <- o33[o33$model == m, ]
  f7 <- d$Rel_ge_0_7_frequency; f39 <- d$Rel_top_0039_frequency
  cat(sprintf("%-15s Rel>=0.7 freq: median %.2f; items >=0.5 freq: %d/%d; items =1.0: %d\n",
    m, median(f7, na.rm = TRUE), sum(f7 >= .5, na.rm = TRUE), nrow(d), sum(f7 >= 1, na.rm = TRUE)))
  cat(sprintf("%-15s top-3.9%% freq: median %.2f; items >=0.5: %d/%d\n",
    m, median(f39, na.rm = TRUE), sum(f39 >= .5, na.rm = TRUE), nrow(d)))
}

cat("\n=== E. selected-subset test quality by coverage (Rel vs alternatives) ===\n")
qc <- qual[qual$coverage_rate %in% c(0.039, 0.087, 0.209), ]
for (m in models) {
  cat(sprintf("\n-- %s --\n", m))
  for (rate in c(0.039, 0.087, 0.209)) {
    d <- qc[qc$model == m & qc$coverage_rate == rate, ]
    cat(sprintf("  coverage %.3f:\n", rate))
    for (mm in c("Rel", "raw_Pearson", "clipped_Pearson", "raw_R2", "training_native_variance")) {
      v <- d$median_test_Pearson[d$method == mm]
      if (!length(v)) next
      cat(sprintf("    %-24s median test Pearson mean %.4f\n", mm, mean(v, na.rm = TRUE)))
    }
  }
}

cat("\n=== F. Jaccard overlap between splits (Rel top-3.9%) ===\n")
jac <- read.csv(file.path(out, "E1_Rel_threshold_pairwise_Jaccard.csv"), check.names = FALSE)
for (rate in unique(jac$coverage_rate)) {
  d <- jac[jac$coverage_rate == rate, ]
  for (m in models) {
    v <- d$Jaccard[d$model == m]
    if (!length(v)) next
    cat(sprintf("coverage %.3f %-15s Jaccard median %.3f range [%.3f, %.3f] (n=%d pairs)\n",
      rate, m, median(v), min(v), max(v), length(v)))
  }
}
