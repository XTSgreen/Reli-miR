# Third-pass: lowest-frequency item and remaining E5 numbers.
out <- "new/results/E1_evaluation"; e5 <- "new/results/E5"
o33 <- read.csv(file.path(out, "E1_original_33_selection_frequency.csv"), check.names = FALSE)
cat("--- lowest original-33 frequency item, by column and model ---\n")
for (m in c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")) {
  d <- o33[o33$model == m, ]
  a <- d[which.min(d$Rel_top_0039_frequency), ]; b <- d[which.min(d$Rel_ge_0_7_frequency), ]
  cat(sprintf("%-15s top-3.9%%: %s (%.2f) | Rel>=0.7: %s (%.2f)\n", m,
    a$miRNA, a$Rel_top_0039_frequency, b$miRNA, b$Rel_ge_0_7_frequency))
}
cat("any frequency exactly 0.10 count:", sum(o33$Rel_top_0039_frequency == 0.10), "\n")
print(o33[o33$Rel_top_0039_frequency <= 0.2 & o33$model == "ElasticNet", c("miRNA", "Rel_top_0039_frequency")])

cat("\n--- E5 ---\n")
md <- read.csv(file.path(e5, "E5_miRNA_identity_expression_score.csv"), check.names = FALSE)
cat("identity_class table:\n"); print(table(md$identity_class))
cat("annotation_status table:\n"); print(table(md$annotation_status))
cat(sprintf("rows with matched UCSC annotation (non-empty ucsc_name): %d\n", sum(nzchar(md$ucsc_name))))
cat(sprintf("rho(Rel, native_expression_mean_366) = %.4f\n",
  cor(md$Rel, md$native_expression_mean_366, method = "spearman", use = "complete.obs")))
cat(sprintf("rho(Rel, native_expression_variance_366) = %.4f\n",
  cor(md$Rel, md$native_expression_variance_366, method = "spearman", use = "complete.obs")))
hi <- md[md$Rel >= 0.7, ]; lo <- md[md$Rel < 0.7, ]
cat(sprintf("high n=%d: median mean %.3f, median var %.3f | other n=%d: median mean %.3f, median var %.3f\n",
  nrow(hi), median(hi$native_expression_mean_366), median(hi$native_expression_variance_366),
  nrow(lo), median(lo$native_expression_mean_366), median(lo$native_expression_variance_366)))
cat(sprintf("high-score with n_host_gene_candidates>0: %d/%d\n",
  sum(hi$n_host_gene_candidates > 0, na.rm = TRUE), nrow(hi)))

hc <- read.csv(file.path(e5, "E5_host_gene_candidates.csv"), check.names = FALSE)
cat(sprintf("host candidate rows %d; distinct miRNA %d; in_model_3000 TRUE rows %d; distinct loci (chrom:start:end:strand) %d\n",
  nrow(hc), length(unique(hc$miRNA)), sum(hc$in_model_3000 %in% c(TRUE, "TRUE")),
  length(unique(paste(hc$chrom, hc$start, hc$end, hc$strand)))))

ha <- read.csv(file.path(e5, "E5_host_expression_associations_all_mRNA.csv"), check.names = FALSE)
cat("all-mRNA association columns:", paste(names(ha), collapse = ", "), "\n")
cat(sprintf("rows %d; distinct miRNA %d\n", nrow(ha), length(unique(ha$miRNA))))
qcol <- grep("^q|padj|FDR|bh", names(ha), ignore.case = TRUE, value = TRUE)
if (length(qcol)) for (q in qcol) cat(sprintf("  %s < 0.05: %d\n", q, sum(ha[[q]] < .05, na.rm = TRUE)))
hm <- read.csv(file.path(e5, "E5_host_expression_associations_model_inputs_only.csv"), check.names = FALSE)
cat(sprintf("model-input rows %d; distinct miRNA %d\n", nrow(hm), length(unique(hm$miRNA))))

cat("\n--- E4 CHOL exact ---\n")
ov <- read.csv("new/results/harmonized/E4/E4_external_overall_metrics.csv", check.names = FALSE)
print(ov, digits = 6)
bs <- read.csv("new/results/harmonized/E4/E4_external_patient_bootstrap_summary.csv", check.names = FALSE)
print(bs[bs$metric %in% c("mean_R2", "median_Pearson", "mean_RMSE", "rho_internal_score_external_Pearson"), ], digits = 6)
th <- read.csv("new/results/harmonized/E4/E4_Rel_threshold_external_coverage.csv", check.names = FALSE)
print(th[th$selection == "Rel_ge_0.7_canonical_33", ], digits = 6)
