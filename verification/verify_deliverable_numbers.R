# Verify every statistic cited in the three deliverable documents against source CSVs.
# Any FAIL must be corrected in the documents before delivery.
out <- "new/results/E1_evaluation"
ok <- 0L; bad <- 0L
chk <- function(label, got, exp, tol = 5e-4) {
  good <- length(got) == length(exp) && all(abs(as.numeric(got) - as.numeric(exp)) <= tol)
  if (isTRUE(good)) { ok <<- ok + 1L; cat(sprintf("PASS  %s\n", label)) }
  else { bad <<- bad + 1L; cat(sprintf("FAIL  %s\n      got: %s\n      exp: %s\n", label,
    paste(round(as.numeric(got), 6), collapse = ", "), paste(round(as.numeric(exp), 6), collapse = ", "))) }
}

perf <- read.csv(file.path(out, "E1_model_performance_by_seed.csv"), check.names = FALSE)
win <- read.csv(file.path(out, "E1_winner_frequency_partial.csv"), check.names = FALSE)
seeds <- sort(unique(perf$seed)); M <- c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")

cat("--- E1 ---\n")
chk("E1 winner counts (EN/PCAEN/PCAR)",
  win$n_seeds[match(M, win$selected_model)], c(14, 6, 0))
for (m in M) {
  d <- perf[perf$model == m, ]
  chk(sprintf("E1 %s median_pearson mean[q05,q95]", m),
    c(mean(d$median_pearson), quantile(d$median_pearson, .05), quantile(d$median_pearson, .95)),
    switch(m, ElasticNet = c(0.47110, 0.44396, 0.49738),
              PCA_ElasticNet = c(0.46433, 0.43875, 0.49214),
              PCA_Ridge = c(0.46529, 0.43746, 0.49464)))
  chk(sprintf("E1 %s mean_r2 mean", m), mean(d$mean_r2),
    switch(m, ElasticNet = 0.20232, PCA_ElasticNet = 0.19637, PCA_Ridge = 0.19775))
  chk(sprintf("E1 %s norm_rmse mean", m), mean(d$normalized_rmse),
    switch(m, ElasticNet = 0.89545, PCA_ElasticNet = 0.89803, PCA_Ridge = 0.89800))
}
d <- perf[perf$model == "Mean", ]
chk("E1 Mean mean_r2 / norm_rmse", c(mean(d$mean_r2), mean(d$normalized_rmse)), c(-0.02692, 1.01943))

pd <- function(a, b, crit = "median_pearson") sapply(seeds, function(s)
  perf[[crit]][perf$seed == s & perf$model == a] - perf[[crit]][perf$seed == s & perf$model == b])
for (p in list(c("ElasticNet", "PCA_ElasticNet"), c("ElasticNet", "PCA_Ridge"), c("PCA_ElasticNet", "PCA_Ridge"))) {
  v <- pd(p[1], p[2])
  chk(sprintf("E1 paired %s-%s mean|n_pos", p[1], p[2]),
    c(mean(v), sum(v > 0)), switch(paste(p, collapse = "-"),
      "ElasticNet-PCA_ElasticNet" = c(0.00677, 17),
      "ElasticNet-PCA_Ridge" = c(0.00581, 16), c(-0.00096, 9)))
}

assoc <- do.call(rbind, lapply(seeds, function(s) read.csv(file.path(out, "by_seed", as.character(s),
  "score_test_associations.csv"), check.names = FALSE)))
for (m in M) chk(sprintf("E1 Rel rho vs RMSE median %s", m),
  median(assoc$spearman_quality_score_vs_test_RMSE_trainingSD[assoc$model == m]),
  switch(m, ElasticNet = -0.57194, PCA_ElasticNet = -0.57323, PCA_Ridge = -0.56042))
chk("E1 Rel rho pooled range",
  range(assoc$spearman_quality_score_vs_test_RMSE_trainingSD[assoc$model == "ElasticNet"]), c(-0.68493, -0.4537))
chk("E1 OOF RMSE rho mean (ElasticNet)",
  mean(assoc$spearman_quality_score_vs_test_RMSE_trainingSD[assoc$model == "ElasticNet" &
    assoc$method == "same_split_training_OOF_RMSE"]), -0.5833)
chk("E1 native mean rho mean (ElasticNet)",
  mean(assoc$spearman_quality_score_vs_test_RMSE_trainingSD[assoc$model == "ElasticNet" &
    assoc$method == "training_native_mean"]), -0.4892)

thr <- do.call(rbind, lapply(seeds, function(s) read.csv(file.path(out, "by_seed", as.character(s),
  "Rel_threshold_counts.csv"), check.names = FALSE)))
thr <- thr[thr$selection_kind == "Rel_threshold", ]
for (m in M) chk(sprintf("E1 median retained at .5/.6/.7 (%s)", m),
  sapply(c(.5, .6, .7), function(t) median(thr$n_selected[thr$model == m & thr$threshold == t])),
  switch(m, ElasticNet = c(199.5, 96, 34), PCA_ElasticNet = c(181, 90.5, 38), PCA_Ridge = c(182.5, 89.5, 37.5)))
chk("E1 Rel>=0.7 retention overall range",
  range(thr$n_selected[thr$threshold == .7]), c(24, 51))

o33 <- read.csv(file.path(out, "E1_original_33_selection_frequency.csv"), check.names = FALSE)
for (m in M) chk(sprintf("E1 orig33 top-3.9%% freq median / n>=0.5 (%s)", m),
  c(median(o33$Rel_top_0039_frequency[o33$model == m]), sum(o33$Rel_top_0039_frequency[o33$model == m] >= .5)),
  c(0.85, if (m == "ElasticNet") 27 else 28))
chk("E1 orig33 Rel>=0.7 freq median EN / in-all-20",
  c(median(o33$Rel_ge_0_7_frequency[o33$model == "ElasticNet"]),
    sum(o33$Rel_ge_0_7_frequency[o33$model == "ElasticNet"] >= 1)), c(0.90, 11))
chk("E1 orig33 min frequency item freq (hsa-mir-490)",
  min(o33$Rel_top_0039_frequency[o33$model == "ElasticNet"]), 0.10)

jac <- read.csv(file.path(out, "E1_Rel_threshold_pairwise_Jaccard.csv"), check.names = FALSE)
chk("E1 Jaccard median range across models/thresholds",
  range(sapply(M, function(m) sapply(c(.5, .6, .7), function(t)
    median(jac$Jaccard[jac$model == m & jac$threshold == t])))), c(0.527, 0.594))

e3 <- do.call(rbind, lapply(seeds, function(s) {
  d <- read.csv(file.path(out, "by_seed", as.character(s), "E3_winner_equal_size_summary.csv"), check.names = FALSE)
  d$seed <- s; d[, c("seed", "subset", "n_miRNA", "median_patient_correlation",
    "correlation_matrix_mean_absolute_difference", "median_variance_ratio",
    "DE_effect_slope_pred_vs_real")] }))
e3$isr <- grepl("^Random_", e3$subset)
dd <- e3[e3$n_miRNA == 33, ]
chk("E3(20) 33 Rel med patient-corr / random mean",
  c(mean(dd$median_patient_correlation[dd$subset == "Rel"]),
    mean(dd$median_patient_correlation[dd$isr])), c(0.8291, 0.4607))
chk("E3(20) 33 Rel CMAD / random mean",
  c(mean(dd$correlation_matrix_mean_absolute_difference[dd$subset == "Rel"]),
    mean(dd$correlation_matrix_mean_absolute_difference[dd$isr])), c(0.1600, 0.2316))
chk("E3(20) 33 Pearson med patient-corr", mean(dd$median_patient_correlation[dd$subset == "Pearson"]), 0.8243)
chk("E3(20) 33 OOF-RMSE med patient-corr",
  mean(dd$median_patient_correlation[dd$subset == "training_OOF_RMSE"]), 0.8347)
for (k in c(33, 74, 177)) {
  dk <- e3[e3$n_miRNA == k & e3$subset == "Rel", ]
  chk(sprintf("E3(20) Rel variance ratio m%d", k), mean(dk$median_variance_ratio, na.rm = TRUE),
    switch(as.character(k), "33" = 0.5142, "74" = 0.4536, 0.3661))
  chk(sprintf("E3(20) Rel DE slope m%d", k), mean(dk$DE_effect_slope_pred_vs_real),
    switch(as.character(k), "33" = 0.5800, "74" = 0.5196, 0.4168))
}

cat("\n--- E2 (fixed split) ---\n")
sp <- read.csv("new/results/E2/E2_complete_score_point_reproduction_check.csv", check.names = FALSE)
cat("E2 score point columns:", paste(names(sp), collapse = ", "), "\n")
pc <- read.csv("new/results/E2/E2_complete_model_paired_contrasts.csv", check.names = FALSE)
cat("E2 paired contrast columns:", paste(names(pc), collapse = ", "), "\n")
print(head(pc, 6))
sb <- read.csv("new/results/E2/E2_complete_sample_bootstrap_summary.csv", check.names = FALSE)
cat("E2 sample summary:\n"); print(sb)

cat("\n--- E3 (fixed split, E3_subset_summary.csv) ---\n")
es <- read.csv("new/results/E3/E3_subset_summary.csv", check.names = FALSE)
row <- function(k, s) es[es$n_miRNA == k & es$subset == s, ][1, ]
a33 <- row(33, "Rel"); a846 <- row(846, "All") 
cat("846 All:", a846$median_patient_correlation, a846$correlation_matrix_correlation,
    a846$median_variance_ratio, a846$Cox_logHR_correlation, "\n")
cat("33 Rel:", a33$median_patient_correlation, a33$correlation_matrix_correlation,
    a33$median_variance_ratio, a33$DE_effect_slope_pred_vs_real, a33$Cox_logHR_correlation,
    a33$correlation_matrix_mean_absolute_difference, "\n")
cat("33 top10/top20:", a33$DE_top10_overlap, a33$DE_top20_overlap, "| top50 =", a33$DE_top50_overlap, "\n")
cat("DE groups:", a33$DE_n_early, a33$DE_n_late, "| Cox:", a33$Cox_events, a33$Cox_censored, "\n")
pci <- read.csv("new/results/E3/E3_paired_comparison_bootstrap_CI.csv", check.names = FALSE)
cat("E3 paired CI columns:", paste(names(pci), collapse = ", "), "\n")
print(head(pci, 8))

cat("\n--- E4 multichort ---\n")
mc <- read.csv("new/results/E4_multicohort/multicohort_E4_summary.csv", check.names = FALSE)
chk("E4multi rho range", range(mc$rho_Pearson_point_sample), c(0.366, 0.587))
chk("E4multi all ranking positive", sum(mc$ranking_transfer_positive), 8)
chk("E4multi all top33 > overall", sum(mc$top33_beats_overall), 8)
chk("E4multi overall median Pearson range", range(mc$overall_median_Pearson), c(0.0647, 0.289))
chk("E4multi top33 median Pearson range", range(mc$top33_median_Pearson), c(0.319, 0.576))

cat("\n--- E5 three-model ---\n")
tm <- read.csv("new/results/E5/E5_threemodel_survival_comparison.csv", check.names = FALSE)
print(tm)
chk("C-M CI", c(tm$CI_low[tm$metric == "C_minus_M"], tm$CI_high[tm$metric == "C_minus_M"]), c(-0.0048, 0.0556))
chk("V-M point", tm$point_estimate[tm$metric == "V_minus_M"], 0.0418)

cat("\n--- shuffled ---\n")
sh <- read.csv("new/results/E1_shuffled/E1_shuffled_20_permutation_summary.csv", check.names = FALSE)
cat("shuffled columns:", paste(names(sh), collapse = ", "), "\n")
print(summary(sh[, grep("mean_r2|rmse", names(sh), ignore.case = TRUE, value = TRUE)]))

cat(sprintf("\n==== checks passed: %d ; failed: %d ====\n", ok, bad))
