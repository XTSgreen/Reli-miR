#!/usr/bin/env Rscript
# Single verification entry point for the Reli-miR open-source release.
#
# Recomputes every statistic quoted in the manuscript response letter directly from the
# tables published under results_tables/, and fails loudly on any mismatch.
#
#   Rscript verification/verify_all.R                 # default tables root: results_tables
#   Rscript verification/verify_all.R path/to/tables  # override
#
# Exit status: 0 when every check passes, 1 when any check fails.
#
# Conventions fixed here (each of these was an ambiguity in earlier ad-hoc checks):
#   * score-versus-error associations are always filtered to method == "Rel";
#   * a model absent from a winner-frequency table counts as zero selections;
#   * "consistently retained" is read from the top-3.9% column, not the Rel>=0.7 column;
#   * host-candidate coverage is reported as 724 miRNAs (distinct genomic loci are 591).

args <- commandArgs(trailingOnly = TRUE)
tables <- if (length(args)) args[1] else "results_tables"
if (!dir.exists(tables)) stop("tables root not found: ", tables)
cat("Verifying against tables root: ", normalizePath(tables, winslash = "/"), "\n\n", sep = "")

n_pass <- 0L; n_fail <- 0L; fails <- character()
# chk: compare numeric vectors elementwise, print the supporting source file on every check.
chk <- function(label, got, exp, src, tol = 5e-4) {
  got <- as.numeric(got); exp <- as.numeric(exp)
  ok <- length(got) == length(exp) && all(is.finite(got) == is.finite(exp)) &&
        all(abs(got - exp) <= tol | (is.na(got) & is.na(exp)))
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
    cat(sprintf("PASS  %s\n", label))
  } else {
    n_fail <<- n_fail + 1L
    fails <<- c(fails, label)
    cat(sprintf("FAIL  %s\n      expected: %s\n      observed: %s\n", label,
                paste(format(exp, digits = 7), collapse = ", "),
                paste(format(got, digits = 7), collapse = ", ")))
  }
  cat(sprintf("      source: %s\n", src))
}
chk_true <- function(label, cond, src) chk(label, if (isTRUE(cond)) 1 else 0, 1, src)
# Row lookup by numeric key with tolerance: stored thresholds can carry
# floating-point residue (e.g. a 0.7 threshold stored as 0.7000000000000001).
pick <- function(d, key, value, tol = 1e-9) {
  i <- which(abs(d[[key]] - value) <= tol)
  if (!length(i)) return(NULL)
  d[i[1], ]
}

p <- function(...) file.path(tables, ...)

## ---------------- E1: 20-split repeated training and evaluation ----------------
cat("\n--- E1 repeated training (20 patient-level splits) ---\n")
win <- read.csv(p("E1", "E1_winner_frequency_partial.csv"), check.names = FALSE)
sel <- function(m) { i <- match(m, win$selected_model); if (is.na(i)) 0 else win$n_seeds[i] }
chk("E1 training-internal selected-model counts (ElasticNet/PCA_ElasticNet/PCA_Ridge)",
    c(sel("ElasticNet"), sel("PCA_ElasticNet"), sel("PCA_Ridge")), c(14, 6, 0),
    p("E1", "E1_winner_frequency_partial.csv"))

perf <- read.csv(p("E1", "E1_model_performance_by_seed.csv"), check.names = FALSE)
seeds <- sort(unique(perf$seed))
exp_dist <- list(ElasticNet = c(0.47110, 0.44396, 0.49738, 0.20232, 0.89545),
                 PCA_ElasticNet = c(0.46433, 0.43875, 0.49214, 0.19637, 0.89803),
                 PCA_Ridge = c(0.46529, 0.43746, 0.49464, 0.19775, 0.89800))
for (m in names(exp_dist)) {
  d <- perf[perf$model == m, ]
  chk(sprintf("E1 %s: median-Pearson mean/q05/q95, mean-R2 mean, norm-RMSE mean", m),
      c(mean(d$median_pearson), quantile(d$median_pearson, .05), quantile(d$median_pearson, .95),
        mean(d$mean_r2), mean(d$normalized_rmse)), exp_dist[[m]],
      p("E1", "E1_model_performance_by_seed.csv"))
}
d <- perf[perf$model == "Mean", ]
chk("E1 Mean baseline: mean R2 and normalised RMSE", c(mean(d$mean_r2), mean(d$normalized_rmse)),
    c(-0.02692, 1.01943), p("E1", "E1_model_performance_by_seed.csv"))

pdiff <- function(a, b) sapply(seeds, function(s)
  perf$median_pearson[perf$seed == s & perf$model == a] -
  perf$median_pearson[perf$seed == s & perf$model == b])
chk("E1 paired median-Pearson difference: ElasticNet - PCA_ElasticNet",
    c(mean(pdiff("ElasticNet", "PCA_ElasticNet")), sum(pdiff("ElasticNet", "PCA_ElasticNet") > 0)),
    c(0.00677, 17), p("E1", "E1_model_performance_by_seed.csv"))
chk("E1 paired median-Pearson difference: ElasticNet - PCA_Ridge",
    c(mean(pdiff("ElasticNet", "PCA_Ridge")), sum(pdiff("ElasticNet", "PCA_Ridge") > 0)),
    c(0.00581, 16), p("E1", "E1_model_performance_by_seed.csv"))
chk("E1 paired median-Pearson difference: PCA_ElasticNet - PCA_Ridge",
    c(mean(pdiff("PCA_ElasticNet", "PCA_Ridge")), sum(pdiff("PCA_ElasticNet", "PCA_Ridge") > 0)),
    c(-0.00096, 9), p("E1", "E1_model_performance_by_seed.csv"))

sa <- read.csv(p("E1", "E1_score_associations_by_seed.csv"), check.names = FALSE)
rel <- function(m, col = "rho_quality_vs_test_RMSE_trainingSD")
  sa[[col]][sa$model == m & sa$method == "Rel"]
chk("E1 reliability (Rel) vs test RMSE: median Spearman per model",
    sapply(c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge"), function(m) median(rel(m))),
    c(-0.57194, -0.57323, -0.56042), p("E1", "E1_score_associations_by_seed.csv"))
chk("E1 reliability (Rel) vs test RMSE: pooled range, ElasticNet",
    range(rel("ElasticNet")), c(-0.68493, -0.45375), p("E1", "E1_score_associations_by_seed.csv"))
chk("E1 reliability (Rel) negative in all 20 splits, per model",
    sapply(c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge"), function(m) sum(rel(m) < 0)),
    c(20, 20, 20), p("E1", "E1_score_associations_by_seed.csv"))
chk("E1 same-split training-OOF RMSE ordering: mean Spearman (ElasticNet)",
    mean(sa$rho_quality_vs_test_RMSE_trainingSD[
      sa$model == "ElasticNet" & sa$method == "same_split_training_OOF_RMSE"]), -0.5833,
    p("E1", "E1_score_associations_by_seed.csv"))
chk("E1 training-set native mean ordering: mean Spearman (ElasticNet)",
    mean(sa$rho_quality_vs_test_RMSE_trainingSD[
      sa$model == "ElasticNet" & sa$method == "training_native_mean"]), -0.4892,
    p("E1", "E1_score_associations_by_seed.csv"))

tc <- read.csv(p("E1", "E1_threshold_counts_by_seed.csv"), check.names = FALSE)
ret <- function(m, t) median(tc$n_selected[tc$model == m & abs(tc$threshold - t) <= 1e-9])
chk("E1 median retained miRNAs at Rel >= 0.5 / 0.6 / 0.7 (per model)",
    c(sapply(c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge"), function(m)
      sapply(c(.5, .6, .7), function(t) ret(m, t)))),
    c(199.5, 96, 34, 181, 90.5, 38, 182.5, 89.5, 37.5),
    p("E1", "E1_threshold_counts_by_seed.csv"))
chk("E1 retained-count range at Rel >= 0.7 (all models, all splits)",
    range(tc$n_selected[abs(tc$threshold - .7) <= 1e-9]), c(24, 51),
    p("E1", "E1_threshold_counts_by_seed.csv"))

o33 <- read.csv(p("E1", "E1_original_33_selection_frequency.csv"), check.names = FALSE)
chk("E1 original-33 frequency in top-3.9% list: median and count >= 0.5 (per model)",
    c(sapply(c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge"), function(m)
      c(median(o33$Rel_top_0039_frequency[o33$model == m]),
        sum(o33$Rel_top_0039_frequency[o33$model == m] >= .5)))),
    c(0.85, 27, 0.85, 28, 0.85, 28), p("E1", "E1_original_33_selection_frequency.csv"))
chk("E1 original-33 retained in ALL 20 splits (top-3.9% column)",
    sapply(c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge"), function(m)
      sum(o33$Rel_top_0039_frequency[o33$model == m] >= 1)), c(11, 10, 10),
    p("E1", "E1_original_33_selection_frequency.csv"))
chk("E1 original-33 frequency in Rel >= 0.7 list: median (ElasticNet)",
    median(o33$Rel_ge_0_7_frequency[o33$model == "ElasticNet"]), 0.90,
    p("E1", "E1_original_33_selection_frequency.csv"))
chk("E1 lowest original-33 retention frequency in any split (ElasticNet)",
    min(o33$Rel_top_0039_frequency[o33$model == "ElasticNet"]), 0.10,
    p("E1", "E1_original_33_selection_frequency.csv"))

js <- read.csv(p("E1", "E1_toplist_jaccard_summary.csv"), check.names = FALSE)
chk("E1 equal-size top-list Jaccard: median range across models and coverages",
    range(js$Jaccard_median), c(0.526, 0.574), p("E1", "E1_toplist_jaccard_summary.csv"), tol = 1e-3)
jth <- read.csv(p("E1", "E1_Rel_threshold_pairwise_Jaccard.csv"), check.names = FALSE)
chk("E1 Rel-threshold list Jaccard: median range across models and thresholds",
    range(sapply(c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge"), function(m)
      sapply(c(.5, .6, .7), function(t) median(jth$Jaccard[jth$model == m & jth$threshold == t]))),
      na.rm = TRUE), c(0.527, 0.594), p("E1", "E1_Rel_threshold_pairwise_Jaccard.csv"), tol = 1e-3)

## ---------------- shuffled permutation control ----------------
cat("\n--- Shuffled pairing control (20 permutations) ---\n")
sh <- read.csv(p("E1_shuffled", "E1_shuffled_20_permutation_summary.csv"), check.names = FALSE)
chk("Shuffled test mean R2 range", range(sh$test_mean_R2), c(-0.02138, -0.01842),
    p("E1_shuffled", "E1_shuffled_20_permutation_summary.csv"))
chk("Shuffled test normalised RMSE range", range(sh$test_normalized_RMSE), c(1.01326, 1.01464),
    p("E1_shuffled", "E1_shuffled_20_permutation_summary.csv"))
chk("Shuffled permutations selecting the first null path (constant prediction)",
    sum(sh$lambda_index == 1), 13, p("E1_shuffled", "E1_shuffled_20_permutation_summary.csv"))

## ---------------- E2: fixed split, scores, thresholds, sample level ----------------
cat("\n--- E2 fixed-split score ablation and calibration ---\n")
sp <- read.csv(p("E2", "E2_complete_score_point_reproduction_check.csv"), check.names = FALSE)
point <- function(s) sp$point_rho_vs_test_RMSE[sp$score == s]
chk("E2 score-versus-test-RMSE Spearman: Rel, raw R2, Pearson",
    c(point("Rel_equal"), point("R2_raw"), point("Pearson")), c(-0.597798, -0.591677, -0.591613),
    p("E2", "E2_complete_score_point_reproduction_check.csv"), tol = 1e-5)
chk("E2 score-versus-test-RMSE Spearman: R2 weights 0.25 / 0.75, native mean, native variance",
    c(point("Rel_R2_weight_025"), point("Rel_R2_weight_075"),
      point("Training_native_mean"), point("Training_native_variance")),
    c(-0.594955, -0.599547, -0.541491, -0.658374),
    p("E2", "E2_complete_score_point_reproduction_check.csv"), tol = 1e-5)

sb <- read.csv(p("E2", "E2_complete_score_joint_bootstrap_summary.csv"), check.names = FALSE)
pe <- sb[sb$score == "Pearson", ]
chk("E2 Pearson minus Rel: point rho difference and 95% interval",
    c(pe$point_rho_baseline_minus_Rel, pe$contrast_CI_q025, pe$contrast_CI_q975),
    c(0.00618446, -0.00703605, 0.02112728),
    p("E2", "E2_complete_score_joint_bootstrap_summary.csv"), tol = 1e-5)

pc <- read.csv(p("E2", "E2_complete_model_paired_contrasts.csv"), check.names = FALSE)
sub <- pc[pc$model_a == "PCA_ElasticNet" & pc$model_b == "PCA_Ridge", ]
val <- function(metric, col) sub[[col]][sub$metric == metric]
chk("E2 PCA_ElasticNet - PCA_Ridge: mean R2 difference and interval",
    c(val("mean_R2", "point_difference"), val("mean_R2", "CI_q025"), val("mean_R2", "CI_q975")),
    c(0.00171388, -0.00256755, 0.00237727),
    p("E2", "E2_complete_model_paired_contrasts.csv"), tol = 1e-6)
chk("E2 PCA_ElasticNet - PCA_Ridge: median Pearson difference and interval",
    c(val("median_Pearson", "point_difference"), val("median_Pearson", "CI_q025"),
      val("median_Pearson", "CI_q975")), c(0.00322380, -0.00192026, 0.00732786),
    p("E2", "E2_complete_model_paired_contrasts.csv"), tol = 1e-6)
chk("E2 PCA_ElasticNet - PCA_Ridge: mean RMSE difference and interval",
    c(val("mean_RMSE", "point_difference"), val("mean_RMSE", "CI_q025"), val("mean_RMSE", "CI_q975")),
    c(-0.00150384, -0.00217734, -0.00058393),
    p("E2", "E2_complete_model_paired_contrasts.csv"), tol = 1e-6)

ss2 <- read.csv(p("E2", "E2_complete_sample_bootstrap_summary.csv"), check.names = FALSE)
sl <- function(s, col) ss2[[col]][ss2$score == s]
chk("E2 sample-level score vs actual test RMSE: rho and interval",
    c(sl("sample_level_reliability", "point_rho_vs_actual_RMSE"),
      sl("sample_level_reliability", "CI_q025"), sl("sample_level_reliability", "CI_q975")),
    c(-0.11926812, -0.3603996, 0.1469063),
    p("E2", "E2_complete_sample_bootstrap_summary.csv"), tol = 1e-6)
chk("E2 sample-level distance component: rho and interval",
    c(sl("distance_score", "point_rho_vs_actual_RMSE"),
      sl("distance_score", "CI_q025"), sl("distance_score", "CI_q975")),
    c(-0.37484726, -0.5762773, -0.1414498),
    p("E2", "E2_complete_sample_bootstrap_summary.csv"), tol = 1e-6)
chk("E2 sample-level model-disagreement component: rho and interval",
    c(sl("model_disagreement", "point_rho_vs_actual_RMSE"),
      sl("model_disagreement", "CI_q025"), sl("model_disagreement", "CI_q975")),
    c(-0.03527558, -0.2829451, 0.2235843),
    p("E2", "E2_complete_sample_bootstrap_summary.csv"), tol = 1e-6)

fn <- read.csv(p("E2", "E2_Rel_threshold_coverage_FNR_point.csv"), check.names = FALSE)
fv <- function(t, col) { r <- pick(fn, "threshold", t); if (is.null(r)) NA_real_ else r[[col]] }
chk("E2 Rel >= 0.7: retained, coverage, precision, recall, predicted FNR",
    c(fv(.7, "n_selected"), fv(.7, "coverage"), fv(.7, "precision"), fv(.7, "recall"), fv(.7, "FNR")),
    c(33, 0.039, 1.00, 0.098, 0.902),
    p("E2", "E2_Rel_threshold_coverage_FNR_point.csv"), tol = 1e-3)
chk("E2 Rel >= 0.6: retained, coverage, precision, recall, predicted FNR",
    c(fv(.6, "n_selected"), fv(.6, "coverage"), fv(.6, "precision"), fv(.6, "recall"), fv(.6, "FNR")),
    c(74, 0.087, 0.986, 0.217, 0.783),
    p("E2", "E2_Rel_threshold_coverage_FNR_point.csv"), tol = 1e-3)
chk("E2 Rel >= 0.5: retained, coverage, precision, recall, predicted FNR",
    c(fv(.5, "n_selected"), fv(.5, "coverage"), fv(.5, "precision"), fv(.5, "recall"), fv(.5, "FNR")),
    c(177, 0.209, 0.842, 0.443, 0.557),
    p("E2", "E2_Rel_threshold_coverage_FNR_point.csv"), tol = 1e-3)

## ---------------- E3: downstream preservation and controls ----------------
cat("\n--- E3 downstream preservation (fixed split and per split) ---\n")
e3 <- read.csv(p("E3", "E3_subset_summary_named.csv"), check.names = FALSE)
rowv <- function(k, s) e3[e3$n_miRNA == k & e3$subset == s, ][1, ]
a <- rowv(846, "All")
chk("E3 all 846 miRNAs: median patient correlation, correlation-matrix preservation, variance ratio, Cox logHR correlation",
    c(a$median_patient_correlation, a$correlation_matrix_correlation, a$median_variance_ratio,
      a$Cox_logHR_correlation), c(0.4250144, 0.5284385, 0.1829140, 0.3802427),
    p("E3", "E3_subset_summary_named.csv"), tol = 1e-6)
b <- rowv(33, "Rel")
chk("E3 Rel>=0.7 (33): median patient correlation, correlation-matrix preservation, variance ratio, Cox logHR correlation",
    c(b$median_patient_correlation, b$correlation_matrix_correlation, b$median_variance_ratio,
      b$Cox_logHR_correlation), c(0.8442343, 0.9520652, 0.5334180, 0.8953377),
    p("E3", "E3_subset_summary_named.csv"), tol = 1e-6)
chk("E3 Rel>=0.7 (33): DE effect slope, correlation-matrix mean absolute difference, all-846 counterpart",
    c(b$DE_effect_slope_pred_vs_real, b$correlation_matrix_mean_absolute_difference,
      a$correlation_matrix_mean_absolute_difference), c(0.6133583, 0.1458938, 0.2410272),
    p("E3", "E3_subset_summary_named.csv"), tol = 1e-6)
chk("E3 Rel>=0.7 (33): top-10 and top-20 DE overlap, top-50 not applicable",
    c(b$DE_top10_overlap, b$DE_top20_overlap, ifelse(is.na(b$DE_top50_overlap), 1, 0)),
    c(0.60, 0.90, 1), p("E3", "E3_subset_summary_named.csv"))
chk("E3 all 846 miRNAs: top-10 and top-20 DE overlap",
    c(a$DE_top10_overlap, a$DE_top20_overlap), c(0.10, 0.15),
    p("E3", "E3_subset_summary_named.csv"))
chk("E3 DE groups and Cox events: early/late patients and events/censored",
    c(b$DE_n_early, b$DE_n_late, b$Cox_events, b$Cox_censored), c(47, 18, 25, 47),
    p("E3", "E3_subset_summary_named.csv"))

ci <- read.csv(p("E3", "E3_paired_comparison_bootstrap_CI.csv"), check.names = FALSE)
g <- function(contrast, size) ci[ci$contrast == contrast & ci$size == size, ]
chk("E3 Rel minus random, patient correlation (33 items): estimate and 95% interval",
    c(g("Rel_minus_Random_patient", 33)$estimate, g("Rel_minus_Random_patient", 33)$lower_95,
      g("Rel_minus_Random_patient", 33)$upper_95), c(0.4007636, 0.3612677, 0.4264983),
    p("E3", "E3_paired_comparison_bootstrap_CI.csv"), tol = 1e-6)
chk("E3 random minus Rel, correlation-matrix mean absolute difference (33 items): estimate and interval",
    c(g("Random_minus_Rel_network", 33)$estimate, g("Random_minus_Rel_network", 33)$lower_95,
      g("Random_minus_Rel_network", 33)$upper_95), c(0.0951621, 0.0831530, 0.1302722),
    p("E3", "E3_paired_comparison_bootstrap_CI.csv"), tol = 1e-6)

ds <- read.csv(p("E3", "E3_downstream_by_split.csv"), check.names = FALSE)
rm_ <- read.csv(p("E3", "E3_random_mean_by_split.csv"), check.names = FALSE)
m33 <- ds[ds$n_miRNA == 33, ]
r33 <- rm_[rm_$n_miRNA == 33, ]
chk("E3 per split, 33 items: Rel and random-mean patient correlation",
    c(mean(m33$median_patient_correlation[m33$subset == "Rel"]),
      mean(r33$random_mean_patient_correlation)),
    c(0.8291, 0.4607), p("E3", "E3_downstream_by_split.csv"), tol = 1e-3)
chk("E3 per split, 33 items: Rel and random-mean correlation-matrix mean absolute difference",
    c(mean(m33$correlation_matrix_mean_absolute_difference[m33$subset == "Rel"]),
      mean(r33$random_mean_correlation_matrix_mean_absolute_difference)),
    c(0.1600, 0.2316), p("E3", "E3_downstream_by_split.csv"), tol = 1e-3)
chk("E3 per split, 33 items: Pearson-subset and training-OOF-subset patient correlation",
    c(mean(m33$median_patient_correlation[m33$subset == "Pearson"]),
      mean(m33$median_patient_correlation[m33$subset == "training_OOF_RMSE"])),
    c(0.8243, 0.8347), p("E3", "E3_downstream_by_split.csv"), tol = 1e-3)
chk("E3 per split: Rel-subset variance ratio and DE effect slope at 33 / 74 / 177 items",
    c(sapply(c(33, 74, 177), function(k) mean(ds$median_variance_ratio[ds$n_miRNA == k & ds$subset == "Rel"])),
      sapply(c(33, 74, 177), function(k) mean(ds$DE_effect_slope_pred_vs_real[ds$n_miRNA == k & ds$subset == "Rel"]))),
    c(0.5142, 0.4536, 0.3661, 0.5800, 0.5196, 0.4168),
    p("E3", "E3_downstream_by_split.csv"), tol = 1e-3)
chk_true("E3 per split: Rel above random-mean patient correlation in every split (33 items)",
    all(m33$median_patient_correlation[m33$subset == "Rel"] > r33$random_mean_patient_correlation),
    p("E3", "E3_random_mean_by_split.csv"))
chk_true("E3 per split: Rel below random-mean correlation-matrix difference in every split (33 items)",
    all(m33$correlation_matrix_mean_absolute_difference[m33$subset == "Rel"] <
        r33$random_mean_correlation_matrix_mean_absolute_difference),
    p("E3", "E3_random_mean_by_split.csv"))

rd <- read.csv(p("E3", "E3_random_subset_distribution_summary.csv"), check.names = FALSE)
chk("E3 fixed-split random-subset distribution: 1,000 lists per size, median patient correlation at 33 items",
    c(rd$n_random_lists[rd$requested_size == 33],
      rd$patient_corr_median[rd$requested_size == 33]), c(1000, 0.44536),
    p("E3", "E3_random_subset_distribution_summary.csv"), tol = 1e-3)
chk("E3 fixed-split random-subset distribution: median correlation-matrix difference at 33 items",
    rd$cmad_median[rd$requested_size == 33], 0.23903,
    p("E3", "E3_random_subset_distribution_summary.csv"), tol = 1e-3)

## ---------------- E4: external transfer, CHOL and 8 cohorts ----------------
cat("\n--- E4 external transfer ---\n")
ov <- read.csv(p("E4_CHOL", "E4_external_overall_metrics.csv"), check.names = FALSE)
chk("E4 CHOL: overlap count, mean R2, median Pearson, mean RMSE",
    c(ov$n_overlap_miRNA, ov$mean_R2, ov$median_Pearson, ov$mean_RMSE),
    c(742, -0.3060348, 0.2892197, 0.8517106),
    p("E4_CHOL", "E4_external_overall_metrics.csv"), tol = 1e-6)
bs <- read.csv(p("E4_CHOL", "E4_complete_external_bootstrap_summary.csv"), check.names = FALSE)
bb <- function(m, col) bs[[col]][bs$metric == m]
chk("E4 CHOL bootstrap: mean R2 interval",
    c(bb("mean_R2", "CI_q025"), bb("mean_R2", "CI_q975")), c(-0.700962, -0.279220),
    p("E4_CHOL", "E4_complete_external_bootstrap_summary.csv"), tol = 1e-5)
chk("E4 CHOL bootstrap: median Pearson interval and valid/requested replicate counts",
    c(bb("median_Pearson", "CI_q025"), bb("median_Pearson", "CI_q975"),
      bb("mean_R2", "B_valid"), bb("mean_R2", "B_requested")),
    c(0.210708, 0.353126, 2000, 2000),
    p("E4_CHOL", "E4_complete_external_bootstrap_summary.csv"), tol = 1e-5)
chk("E4 CHOL point estimates: mean R2, median Pearson, mean RMSE",
    c(bb("mean_R2", "point_estimate"), bb("median_Pearson", "point_estimate"),
      bb("mean_RMSE", "point_estimate")), c(-0.306035, 0.289220, 0.851711),
    p("E4_CHOL", "E4_complete_external_bootstrap_summary.csv"), tol = 1e-5)
chk("E4 CHOL bootstrap: Spearman of internal Rel vs external Pearson, point estimate and interval",
    c(bb("rho_internal_Rel_vs_Pearson", "point_estimate"),
      bb("rho_internal_Rel_vs_Pearson", "CI_q025"), bb("rho_internal_Rel_vs_Pearson", "CI_q975")),
    c(0.366358, 0.224371, 0.381486),
    p("E4_CHOL", "E4_complete_external_bootstrap_summary.csv"), tol = 1e-5)
th <- read.csv(p("E4_CHOL", "E4_Rel_threshold_external_coverage.csv"), check.names = FALSE)
t33 <- th[th$selection == "Rel_ge_0.7_canonical_33", ]
chk("E4 CHOL Rel>=0.7 list: external overlap, median R2, median Pearson, reconstructable count",
    c(t33$external_overlap, t33$median_R2, t33$median_Pearson, t33$reconstructable_n),
    c(32, 0.0494735, 0.576221, 16),
    p("E4_CHOL", "E4_Rel_threshold_external_coverage.csv"), tol = 1e-5)

mc <- read.csv(p("E4_multicohort", "multicohort_E4_summary.csv"), check.names = FALSE)
chk("E4 multicohort: number of cohorts and total paired patients",
    c(nrow(mc), sum(mc$n_patients)), c(8, 3269),
    p("E4_multicohort", "multicohort_E4_summary.csv"))
chk("E4 multicohort: sample point-estimate rho range for Rel vs external Pearson",
    range(mc$rho_Pearson_point_sample), c(0.366, 0.587),
    p("E4_multicohort", "multicohort_E4_summary.csv"), tol = 1e-3)
chk("E4 multicohort: cohorts with positive ranking transfer and with top-33 above overall",
    c(sum(mc$ranking_transfer_positive), sum(mc$top33_beats_overall)), c(8, 8),
    p("E4_multicohort", "multicohort_E4_summary.csv"))
chk_true("E4 multicohort: lower 95% bound of the ranking-transfer rho above zero in all cohorts",
    all(mc$rho_Pearson_low > 0), p("E4_multicohort", "multicohort_E4_summary.csv"))
chk("E4 multicohort: overall and top-33 median Pearson ranges",
    c(range(mc$overall_median_Pearson), range(mc$top33_median_Pearson)),
    c(0.0647, 0.289, 0.319, 0.576),
    p("E4_multicohort", "multicohort_E4_summary.csv"), tol = 1e-3)

## ---------------- E5: identity, host analysis, three-model comparison ----------------
cat("\n--- E5 biological annotation and prognostic comparison ---\n")
tm <- read.csv(p("E5", "E5_threemodel_survival_comparison.csv"), check.names = FALSE)
tv <- function(m, col) tm[[col]][tm$metric == m]
chk("E5 three-model test C-index: mRNA, virtual miRNA, combined, real-miRNA diagnostic",
    c(tv("M_mRNA", "point_estimate"), tv("V_virtual", "point_estimate"),
      tv("C_combined", "point_estimate"), tv("R_real_diag", "point_estimate")),
    c(0.6190476, 0.6608358, 0.6423712, 0.6482021),
    p("E5", "E5_threemodel_survival_comparison.csv"), tol = 1e-6)
chk("E5 combined minus mRNA C-index difference and 95% interval (includes zero)",
    c(tv("C_minus_M", "point_estimate"), tv("C_minus_M", "CI_low"), tv("C_minus_M", "CI_high")),
    c(0.0233236, -0.0048311, 0.0556179),
    p("E5", "E5_threemodel_survival_comparison.csv"), tol = 1e-6)
chk("E5 virtual minus mRNA C-index difference and 95% interval",
    c(tv("V_minus_M", "point_estimate"), tv("V_minus_M", "CI_low"), tv("V_minus_M", "CI_high")),
    c(0.0417881, -0.0244456, 0.1049412),
    p("E5", "E5_threemodel_survival_comparison.csv"), tol = 1e-6)

md <- read.csv(p("E5", "E5_miRNA_identity_expression_score.csv"), check.names = FALSE)
hi <- md[md$Rel >= 0.7, ]; lo <- md[md$Rel < 0.7, ]
chk("E5 identity table: total miRNAs and high-reliability miRNAs",
    c(nrow(md), nrow(hi)), c(846, 33), p("E5", "E5_miRNA_identity_expression_score.csv"))
chk("E5 Spearman(Rel, expression mean) and Spearman(Rel, expression variance)",
    c(cor(md$Rel, md$native_expression_mean_366, method = "spearman", use = "complete.obs"),
      cor(md$Rel, md$native_expression_variance_366, method = "spearman", use = "complete.obs")),
    c(0.6942, 0.7617), p("E5", "E5_miRNA_identity_expression_score.csv"), tol = 1e-3)
chk("E5 median expression mean/variance for high-reliability and remaining miRNAs",
    c(median(hi$native_expression_mean_366), median(hi$native_expression_variance_366),
      median(lo$native_expression_mean_366), median(lo$native_expression_variance_366)),
    c(7.408, 2.052, 0.702, 0.405), p("E5", "E5_miRNA_identity_expression_score.csv"), tol = 1e-3)
chk("E5 high-reliability miRNAs with at least one candidate host gene",
    sum(hi$n_host_gene_candidates > 0, na.rm = TRUE), 28,
    p("E5", "E5_miRNA_identity_expression_score.csv"))

rp <- read.csv(p("E5", "E5_host_expression_stratified_label_permutation.csv"), check.names = FALSE)
allv <- rp[rp$analysis_scope == "all_18204_mRNA", ]
chk("E5 host co-expression permutation (all 18,204 mRNAs): strata, common support, weighted difference, p",
    c(allv$n_common_support_strata, allv$n_high_common_support, allv$n_background_common_support,
      allv$observed_stratum_weighted_mean_abs_rho_difference, allv$label_permutation_p_two_sided),
    c(8, 13, 40, 0.0270870, 0.647135),
    p("E5", "E5_host_expression_stratified_label_permutation.csv"), tol = 1e-5)
modv <- rp[rp$analysis_scope == "model_inputs_only_3000", ]
chk("E5 host co-expression permutation (model-input genes only): common support, weighted difference, p",
    c(modv$n_high_common_support, modv$n_background_common_support,
      modv$observed_stratum_weighted_mean_abs_rho_difference, modv$label_permutation_p_two_sided),
    c(2, 6, -0.0465422, 0.473153),
    p("E5", "E5_host_expression_stratified_label_permutation.csv"), tol = 1e-5)

## ---------------- verdict ----------------
cat(sprintf("\n==== checks passed: %d ; failed: %d ====\n", n_pass, n_fail))
if (n_fail > 0L) {
  cat("Failed checks:\n"); for (f in fails) cat("  -", f, "\n")
  cat("Note: a FAIL means a published number does not reproduce from the published tables.\n")
  quit(status = 1L)
}
cat("All published statistics reproduce from results_tables/.\n")
quit(status = 0L)
