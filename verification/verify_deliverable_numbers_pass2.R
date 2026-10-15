# Second-pass verification: exact values for the disputed claims.
out <- "new/results/E1_evaluation"; seeds <- 20261014:20261033
M <- c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")

cat("--- (a) Rel-only rho vs test RMSE ---\n")
assoc <- do.call(rbind, lapply(seeds, function(s) read.csv(file.path(out, "by_seed", as.character(s),
  "score_test_associations.csv"), check.names = FALSE)))
for (m in M) {
  v <- assoc$spearman_quality_score_vs_test_RMSE_trainingSD[assoc$model == m & assoc$method == "Rel"]
  cat(sprintf("%-15s Rel-only: mean %+.5f median %+.5f range [%+.5f, %+.5f]; negative in %d/%d\n",
    m, mean(v), median(v), min(v), max(v), sum(v < 0), length(v)))
}

cat("\n--- (b) original-33 items retained in ALL 20 splits ---\n")
o33 <- read.csv(file.path(out, "E1_original_33_selection_frequency.csv"), check.names = FALSE)
for (m in M) {
  d <- o33[o33$model == m, ]
  cat(sprintf("%-15s top-3.9%% freq==1: %d/33 (>=0.5: %d) | Rel>=0.7 freq==1: %d/33 (>=0.5: %d)\n",
    m, sum(d$Rel_top_0039_frequency >= 1), sum(d$Rel_top_0039_frequency >= .5),
    sum(d$Rel_ge_0_7_frequency >= 1), sum(d$Rel_ge_0_7_frequency >= .5)))
  lo <- d[which.min(d$Rel_top_0039_frequency), ]
  cat(sprintf("   lowest item: %s freq %.2f\n", lo$miRNA, lo$Rel_top_0039_frequency))
}

cat("\n--- (c) E2 fixed-split paired contrasts (PCA_ElasticNet vs PCA_Ridge) ---\n")
pc <- read.csv("new/results/E2/E2_complete_model_paired_contrasts.csv", check.names = FALSE)
sub <- pc[pc$model_a == "PCA_ElasticNet" & pc$model_b == "PCA_Ridge", c("metric", "point_difference", "CI_q025", "CI_q975")]
print(sub, digits = 6)

cat("\n--- (d) shuffled: exact ranges and constant-prediction count ---\n")
sh <- read.csv("new/results/E1_shuffled/E1_shuffled_20_permutation_summary.csv", check.names = FALSE)
cat(sprintf("test_mean_R2 range: %.5f to %.5f\n", min(sh$test_mean_R2), max(sh$test_mean_R2)))
cat(sprintf("test_normalized_RMSE range: %.5f to %.5f\n", min(sh$test_normalized_RMSE), max(sh$test_normalized_RMSE)))
cat(sprintf("lambda_index==1 (first null path) count: %d/20\n", sum(sh$lambda_index == 1)))
cat(sprintf("test_median_pearson non-finite count: %d/20\n", sum(!is.finite(sh$test_median_pearson))))

cat("\n--- (e) E5 numbers ---\n")
e5 <- "new/results/E5"
mp <- read.csv(file.path(e5, "E5_locus_gene_identity_all_overlaps.csv"), check.names = FALSE)
cat("locus identity rows:", nrow(mp), "| columns:", paste(names(mp), collapse = ", "), "\n")
hc <- read.csv(file.path(e5, "E5_host_gene_candidates.csv"), check.names = FALSE)
cat("host candidates rows:", nrow(hc), "| columns:", paste(names(hc), collapse = ", "), "\n")
md <- read.csv(file.path(e5, "E5_miRNA_identity_expression_score.csv"), check.names = FALSE)
cat("identity/expression rows:", nrow(md), "| columns:", paste(names(md), collapse = ", "), "\n")
cat("high-score (Rel>=0.7) rows:", sum(md$Rel >= 0.7, na.rm = TRUE), "\n")
for (cc in c("mean", "variance", "var", "logvar")) if (cc %in% names(md)) {
  v <- md[[cc]]; r <- md$Rel
  cat(sprintf("cor(Rel, %s) = %.4f\n", cc, cor(r, v, method = "spearman", use = "complete.obs")))
}
rp <- read.csv(file.path(e5, "E5_host_expression_stratified_label_permutation.csv"), check.names = FALSE)
cat("permutation file columns:", paste(names(rp), collapse = ", "), "\n"); print(rp, digits = 6)
ha <- read.csv(file.path(e5, "E5_host_expression_associations_all_mRNA.csv"), check.names = FALSE)
cat("host associations rows:", nrow(ha), "| BH<0.05:", if ("q_value" %in% names(ha)) sum(ha$q_value < .05, na.rm = TRUE) else
  if ("qvalue" %in% names(ha)) sum(ha$qvalue < .05, na.rm = TRUE) else "n/a", "\n")
hmi <- read.csv(file.path(e5, "E5_host_expression_associations_model_inputs_only.csv"), check.names = FALSE)
cat("model-input host associations rows:", nrow(hmi), "\n")
