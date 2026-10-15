# Descriptive exports of the audited E0 input; no model training or selection.
out <- "new/results/E0"
can <- readRDS(file.path(out, "canonical_input.rds"))
stopifnot(length(can$gene_names) == 3000L, length(can$miRNA_names) == 846L)
write.csv(data.frame(gene = can$gene_names,
  training_center = can$x_center[can$gene_names], training_SD = can$x_scale[can$gene_names],
  training_imputation_median = can$x_median[can$gene_names],
  source = "audited original PCA_ElasticNet model input; training-side variance selection"),
  file.path(out, "E0_model_input_gene_list.csv"), row.names = FALSE)
Y <- rbind(can$Y_train_native, can$Y_cal_native, can$Y_test_native)
stopifnot(nrow(Y) == 366L, !anyDuplicated(rownames(Y)))
write.csv(data.frame(miRNA = colnames(Y), n_finite = colSums(is.finite(Y)),
  input_mean_366 = colMeans(Y), input_SD_366 = apply(Y, 2, sd),
  input_median_366 = apply(Y, 2, median), input_min_366 = apply(Y, 2, min), input_max_366 = apply(Y, 2, max),
  training_center = can$y_center[colnames(Y)], training_SD = can$y_scale[colnames(Y)],
  training_imputation_median = can$y_median[colnames(Y)],
  scale_note = "provided input expression scale; physical units and original normalization undocumented"),
  file.path(out, "E0_miRNA_input_descriptive_statistics.csv"), row.names = FALSE)
clinical <- rbind(can$clinical_train, can$clinical_cal, can$clinical_test)
stopifnot(nrow(clinical) == 366L, !anyDuplicated(clinical$patient_id))
rows <- list()
for (part in c("all", "train", "cal", "test")) {
  cc <- if (part == "all") clinical else can[[paste0("clinical_", part)]]
  for (v in c("OS_status", "stage", "gender")) {
    raw <- as.character(cc[[v]]); raw[is.na(raw) | !nzchar(raw)] <- "Missing"
    t <- table(raw)
    rows[[length(rows) + 1L]] <- data.frame(split = part, variable = v, level = names(t), n = as.integer(t), split_n = nrow(cc))
  }
}
write.csv(do.call(rbind, rows), file.path(out, "E0_clinical_group_counts.csv"), row.names = FALSE)
write.csv(data.frame(patient_id = clinical$patient_id,
  split = rep(c("train", "calibration", "test"), c(222L, 72L, 72L))),
  file.path(out, "E0_fixed_split_patient_manifest.csv"), row.names = FALSE)
nn <- readRDS("output/models/MLP_model.rds")
stopifnot(nn$type == "mlp", nn$npc == 50, nn$size == 30)
write.csv(data.frame(parameter = c("implementation", "saved_PCA_inputs", "saved_hidden_units", "saved_network_n", "linear_output", "skip_connections", "weight_decay", "max_iterations", "maximum_weights", "randomness", "batch_size_learning_rate_early_stopping"),
  value = c("R nnet::nnet; code/StrongMiRRecon_v2.R", nn$npc, nn$size, paste(nn$fit$n, collapse = "/"), "TRUE", "TRUE", "1e-4", "300", "500000", "entry set.seed(123), shared RNG stream; per-model initialization seed not separately recorded", "not specified in nnet call; do not invent PyTorch training parameters"),
  provenance = c("saved model and source call", rep("saved model", 3), rep("source call, code/StrongMiRRecon_v2.R:1050", 5), "script entry", "source call")),
  file.path(out, "E0_original_MLP_training_settings.csv"), row.names = FALSE)
cat("Data description passed: 3000 genes, 846 miRNAs, 366 clinical records; saved MLP 50/30/846.\n")
