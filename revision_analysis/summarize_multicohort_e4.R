# Aggregate the 8-cohort E4 transfer results (CHOL + 7 TCGA cohorts) into
# one summary table and a direction-consistency check.
suppressWarnings(suppressMessages(library(readr)))
cohorts <- c("CHOL", "KIRC", "LUAD", "COAD", "STAD", "PRAD", "THCA", "HNSC")
root <- function(co) if (co == "CHOL") "new/results/harmonized" else
  file.path("new/results/harmonized_multicohort", co)
rows <- lapply(cohorts, function(co) {
  d <- root(co)
  ov <- read_csv(file.path(d, "E4", "E4_external_overall_metrics.csv"), show_col_types = FALSE)
  meta <- read_csv(file.path(d, "E4", "E4_run_metadata.csv"), show_col_types = FALSE)
  n_pat <- if ("n_patients" %in% names(ov)) ov$n_patients[1] else meta$n_patients[1]
  thr <- read_csv(file.path(d, "E4", "E4_Rel_threshold_external_coverage.csv"), show_col_types = FALSE)
  bs <- read_csv(file.path(d, "E4", "E4_external_patient_bootstrap_summary.csv"), show_col_types = FALSE)
  t33 <- thr[thr$selection == "Rel_ge_0.7_canonical_33", ]
  rhoP <- bs[bs$metric == "rho_internal_score_external_Pearson", ]
  rhoR <- bs[bs$metric == "rho_internal_score_external_RMSE", ]
  rhoR2 <- bs[bs$metric == "rho_internal_score_external_R2", ]
  by_mir <- read_csv(file.path(d, "E4", "E4_external_metrics_by_miRNA.csv"), show_col_types = FALSE)
  ok <- is.finite(by_mir$internal_Rel) & is.finite(by_mir$Pearson)
  rhoP_point <- suppressWarnings(cor(by_mir$internal_Rel[ok], by_mir$Pearson[ok], method = "spearman"))
  data.frame(cohort = co, n_patients = n_pat,
    n_overlap_miRNA = ov$n_overlap_miRNA[1],
    overall_mean_R2 = ov$mean_R2[1], overall_median_Pearson = ov$median_Pearson[1],
    top33_median_Pearson = t33$median_Pearson[1], top33_median_R2 = t33$median_R2[1],
    top33_reconstructable_n = t33$reconstructable_n[1], top33_external_n = t33$external_overlap[1],
    rho_Pearson_point_sample = rhoP_point,
    rho_Pearson_point = rhoP$estimate_median[1], rho_Pearson_low = rhoP$CI_low[1], rho_Pearson_high = rhoP$CI_high[1],
    rho_RMSE_point = rhoR$estimate_median[1], rho_RMSE_low = rhoR$CI_low[1], rho_RMSE_high = rhoR$CI_high[1],
    rho_R2_point = rhoR2$estimate_median[1], rho_R2_low = rhoR2$CI_low[1], rho_R2_high = rhoR2$CI_high[1],
    stringsAsFactors = FALSE)
})
s <- do.call(rbind, rows)
s$ranking_transfer_positive <- s$rho_Pearson_low > 0
s$top33_beats_overall <- s$top33_median_Pearson > s$overall_median_Pearson
write.csv(s, "new/results/E4_multicohort/multicohort_E4_summary.csv", row.names = FALSE)
print(s[, c("cohort", "n_patients", "overall_mean_R2", "overall_median_Pearson",
            "top33_median_Pearson", "rho_Pearson_point", "rho_Pearson_low",
            "rho_Pearson_high", "ranking_transfer_positive", "top33_beats_overall")])
cat("cohorts with positive ranking transfer:", sum(s$ranking_transfer_positive), "/", nrow(s), "\n")
cat("cohorts where top-33 beats overall median Pearson:", sum(s$top33_beats_overall), "/", nrow(s), "\n")
