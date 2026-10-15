#!/usr/bin/env Rscript
# E2 scoring and E4 external-transfer analyses for the canonical PCA_ElasticNet run.
# The script deliberately requires the E0 canonical RDS; it never silently falls
# back to older summary files or changes scores using test/external outcomes.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
write_csv_safe <- function(x, path) readr::write_csv(as.data.frame(x), path, na = "")

assert_matrix <- function(x, label) {
  if (is.null(x) || length(dim(x)) != 2L) stop("canonical input missing matrix: ", label)
  x <- as.matrix(x); storage.mode(x) <- "numeric"; x
}

metric_cols <- function(real, pred) {
  stopifnot(identical(dim(real), dim(pred)))
  n <- ncol(real)
  r2 <- pearson <- rmse <- rep(NA_real_, n)
  for (j in seq_len(n)) {
    ok <- is.finite(real[, j]) & is.finite(pred[, j])
    a <- real[ok, j]; b <- pred[ok, j]
    if (length(a) >= 3L) {
      sst <- sum((a - mean(a))^2)
      if (sst > 0) r2[j] <- 1 - sum((a - b)^2) / sst
      if (sd(a) > 0 && sd(b) > 0) pearson[j] <- suppressWarnings(cor(a, b))
      rmse[j] <- sqrt(mean((a - b)^2))
    }
  }
  tibble(miRNA = colnames(real), R2 = r2, Pearson = pearson, RMSE = rmse)
}

rank_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

rank_score <- function(sc, names) {
  # Stable tie-breaking by miRNA name, never by evaluation outcomes.
  names[order(-sc[names], names, na.last = NA)]
}

score_variants <- function(cal, train_y, train_native = train_y) {
  nm <- cal$miRNA
  r2 <- setNames(cal$R2, nm); pr <- setNames(cal$Pearson, nm)
  rm <- setNames(cal$RMSE, nm)
  sy <- apply(train_native, 2, sd, na.rm = TRUE)
  my <- colMeans(train_native, na.rm = TRUE)
  list(
    Rel_equal = pmin(pmax(r2, 0), 1) * 0.5 + pmin(pmax(pr, 0), 1) * 0.5,
    R2_raw = r2,
    R2_clipped = pmin(pmax(r2, 0), 1),
    Pearson = pr,
    Pearson_clipped = pmin(pmax(pr, 0), 1),
    Rel_R2_weight_025 = pmin(pmax(r2, 0), 1) * 0.25 + pmin(pmax(pr, 0), 1) * 0.75,
    Rel_R2_weight_075 = pmin(pmax(r2, 0), 1) * 0.75 + pmin(pmax(pr, 0), 1) * 0.25,
    Calibration_RMSE_rank = -rm,
    Training_native_mean = my,
    Training_native_variance = sy^2
  )
}

bootstrap_assoc <- function(cal_y, cal_p, test_y, test_p, train_native,
                            B = 2000L, seed = 20261014L) {
  set.seed(seed)
  nm <- colnames(cal_y)
  out <- vector("list", B)
  for (b in seq_len(B)) {
    itr <- sample.int(nrow(train_native), nrow(train_native), replace = TRUE)
    ic <- sample.int(nrow(cal_y), nrow(cal_y), replace = TRUE)
    it <- sample.int(nrow(test_y), nrow(test_y), replace = TRUE)
    cm <- metric_cols(cal_y[ic, , drop = FALSE], cal_p[ic, , drop = FALSE])
    tm <- metric_cols(test_y[it, , drop = FALSE], test_p[it, , drop = FALSE])
    vars <- score_variants(cm, cal_y[ic, , drop = FALSE], train_native[itr, , drop = FALSE])
    out[[b]] <- tibble(
      replicate = b,
      Rel_equal_vs_test_RMSE = rank_cor(vars$Rel_equal[nm], tm$RMSE),
      Pearson_vs_test_RMSE = rank_cor(vars$Pearson[nm], tm$RMSE),
      R2_vs_test_RMSE = rank_cor(vars$R2_clipped[nm], tm$RMSE),
      Rel_equal_vs_test_R2 = rank_cor(vars$Rel_equal[nm], tm$R2),
      Rel_equal_vs_test_Pearson = rank_cor(vars$Rel_equal[nm], tm$Pearson)
    )
  }
  bind_rows(out)
}

threshold_metrics <- function(scores, test_metrics, thresholds, label_rule) {
  nm <- intersect(names(scores), test_metrics$miRNA)
  t <- test_metrics[match(nm, test_metrics$miRNA), ]
  y <- label_rule(t$R2, t$Pearson)
  bind_rows(lapply(thresholds, function(cut) {
    keep <- is.finite(scores[nm]) & scores[nm] >= cut
    tp <- sum(keep & y, na.rm = TRUE); fp <- sum(keep & !y, na.rm = TRUE)
    fn <- sum(!keep & y, na.rm = TRUE)
    tibble(threshold = cut, n_total = length(nm), n_selected = sum(keep),
           coverage = mean(keep), n_reconstructable = sum(y, na.rm = TRUE),
           precision = if (tp + fp) tp / (tp + fp) else NA_real_,
           recall = if (tp + fn) tp / (tp + fn) else NA_real_,
           FNR = if (tp + fn) fn / (tp + fn) else NA_real_,
           valid_denominator = tp + fn)
  }))
}

bootstrap_fnr <- function(cal_y, cal_p, test_y, test_p, B, score_thresholds,
                          r_thresholds = c(.3, .5, .7), seed = 20261015L) {
  set.seed(seed); out <- vector("list", B)
  for (b in seq_len(B)) {
    ic <- sample.int(nrow(cal_y), nrow(cal_y), replace = TRUE)
    it <- sample.int(nrow(test_y), nrow(test_y), replace = TRUE)
    cm <- metric_cols(cal_y[ic, , drop = FALSE], cal_p[ic, , drop = FALSE])
    tm <- metric_cols(test_y[it, , drop = FALSE], test_p[it, , drop = FALSE])
    sc <- score_variants(cm, cal_y[ic, , drop = FALSE])$Rel_equal
    nm <- intersect(names(sc), tm$miRNA); tm <- tm[match(nm, tm$miRNA), ]
    for (r_threshold in r_thresholds) for (cut in score_thresholds) {
      pos <- tm$R2 > 0 & tm$Pearson >= r_threshold
      keep <- is.finite(sc[nm]) & sc[nm] >= cut
      tp <- sum(keep & pos, na.rm = TRUE); fp <- sum(keep & !pos, na.rm = TRUE)
      fn <- sum(!keep & pos, na.rm = TRUE)
      out[[b]] <- rbind(out[[b]], data.frame(replicate = b, threshold = cut,
        Pearson_label_threshold = r_threshold,
        n_positive = sum(pos, na.rm = TRUE), valid_denominator = tp + fn,
        precision = if (tp + fp) tp/(tp + fp) else NA_real_,
        recall = if (tp + fn) tp/(tp + fn) else NA_real_,
        FNR = if (tp + fn) fn/(tp + fn) else NA_real_))
    }
  }
  bind_rows(out)
}

summarise_boot <- function(x, keys) {
  x %>% group_by(across(all_of(keys))) %>% summarise(
    B_valid_FNR = sum(is.finite(FNR)), B_valid_precision = sum(is.finite(precision)),
    B_valid_recall = sum(is.finite(recall)),
    positive_median = median(n_positive, na.rm = TRUE),
    precision_lo = quantile(precision, .025, na.rm = TRUE, names = FALSE),
    precision_hi = quantile(precision, .975, na.rm = TRUE, names = FALSE),
    recall_lo = quantile(recall, .025, na.rm = TRUE, names = FALSE),
    recall_hi = quantile(recall, .975, na.rm = TRUE, names = FALSE),
    FNR_lo = quantile(FNR, .025, na.rm = TRUE, names = FALSE),
    FNR_hi = quantile(FNR, .975, na.rm = TRUE, names = FALSE),
    FNR_NA_fraction = mean(!is.finite(FNR)), .groups = "drop")
}

sample_level_analysis <- function(x, output_dir) {
  # Recompute both sample-level components from the canonical split and the
  # registered, same-split model predictions. Preserve 19/20 only as a separate
  # provenance comparison, never as the primary reconstructed score.
  obj <- x$metadata_lightweight_object
  if (is.null(obj) && !is.null(x$lightweight_object_path) && file.exists(x$lightweight_object_path))
    obj <- readRDS(x$lightweight_object_path)
  if (is.null(obj)) {
    write_csv_safe(tibble(status = "pending", reason = "canonical metadata_lightweight_object absent"),
                   file.path(output_dir, "E2_sample_level_status.csv")); return(invisible(NULL))
  }
  pca_x <- x$metadata_pca_distance %||% obj$pca_distance
  if (is.null(pca_x) || is.null(x$X_train_z)) {
    write_csv_safe(tibble(status = "pending", reason = "distance PCA or canonical X_train_z absent"),
                   file.path(output_dir, "E2_sample_level_status.csv")); return(invisible(NULL))
  }
  n_pc <- min(obj$n_pc_distance %||% 20L, ncol(x$X_train_z), ncol(pca_x$rotation))
  train_pc <- predict(pca_x, newdata = x$X_train_z)[, seq_len(n_pc), drop = FALSE]
  center <- colMeans(train_pc)
  variance <- apply(train_pc, 2, var)
  variance[!is.finite(variance) | variance < 1e-8] <- 1
  d_train <- rowSums(sweep(sweep(train_pc, 2, center, "-"), 2, sqrt(variance), "/")^2)
  distance_q95 <- as.numeric(quantile(d_train, .95, na.rm = TRUE))
  if (!is.finite(distance_q95) || distance_q95 <= 0) distance_q95 <- median(d_train, na.rm = TRUE)
  distfun <- function(xx) {
    pc <- predict(pca_x, newdata = xx)[, seq_len(n_pc), drop = FALSE]
    d <- rowSums(sweep(sweep(pc, 2, center, "-"), 2, sqrt(variance), "/")^2)
    exp(-d / distance_q95)
  }
  err <- sqrt(rowMeans((x$Y_test_z - x$pred_test_z)^2, na.rm = TRUE))
  ids <- rownames(x$Y_test_z)
  candidate <- c("Ridge", "ElasticNet", "PLS", "PCA_Ridge", "PCA_ElasticNet", "MLP", "Pathway_Ridge")
  reg <- obj$prediction_registry
  cal_predictions <- test_predictions <- list()
  if (is.null(reg)) {
    write_csv_safe(tibble(status = "pending", reason = "prediction registry absent; model disagreement cannot be reconstructed"),
                   file.path(output_dir, "E2_sample_level_status.csv")); return(invisible(NULL))
  }
  for (m in intersect(candidate, reg$model)) {
    row <- reg[match(m, reg$model), ]
    cf <- as.character(row$validation_prediction_file)
    tf <- as.character(row$test_prediction_file)
    if (file.exists(cf) && file.exists(tf)) {
      cal_predictions[[m]] <- readRDS(cf)
      test_predictions[[m]] <- readRDS(tf)
    }
  }
  if (length(cal_predictions) < 2L || length(test_predictions) < 2L) {
    write_csv_safe(tibble(status = "pending", reason = "fewer than two registered model prediction pairs available"),
                   file.path(output_dir, "E2_sample_level_status.csv")); return(invisible(NULL))
  }
  align_pred <- function(preds, target_y) {
    ids0 <- rownames(target_y); vars <- colnames(target_y)
    lapply(preds, function(p) {
      p <- as.matrix(p)
      if (!all(ids0 %in% rownames(p)) || !all(vars %in% colnames(p))) stop("registered predictions do not align with canonical IDs/features")
      p[ids0, vars, drop = FALSE]
    })
  }
  cal_predictions <- align_pred(cal_predictions, x$Y_cal_z)
  test_predictions <- align_pred(test_predictions, x$Y_test_z)
  if ("PCA_ElasticNet" %in% names(test_predictions) &&
      max(abs(test_predictions$PCA_ElasticNet - x$pred_test_z), na.rm = TRUE) > 1e-8) {
    stop("registered PCA_ElasticNet test prediction does not match canonical pred_test_z")
  }
  disagreement <- function(preds) {
    ar <- simplify2array(preds)
    apply(ar, 1, function(z) mean(apply(z, 1, sd, na.rm = TRUE), na.rm = TRUE))
  }
  dcal <- disagreement(cal_predictions); dtest <- disagreement(test_predictions)
  dmin <- min(dcal, na.rm = TRUE); dmax <- max(dcal, na.rm = TRUE)
  agreement <- if (is.finite(dmin) && is.finite(dmax) && dmax > dmin)
    pmin(pmax(1 - (dtest - dmin)/(dmax - dmin), 0), 1) else rep(NA_real_, length(dtest))
  test_dist <- distfun(x$X_test_z)
  rel <- rowMeans(cbind(test_dist, agreement), na.rm = TRUE)
  legacy <- obj$sample_reliability_test
  dat <- tibble(patient_id = ids, dataset = "test",
    distance_score = test_dist, model_disagreement = dtest, agreement_score = agreement,
    sample_level_reliability = rel, PCA_ElasticNet_sample_RMSE = err,
    legacy_saved_reliability = legacy$sample_level_reliability[match(ids, legacy$patient_id)],
    legacy_saved_disagreement = legacy$model_disagreement[match(ids, legacy$patient_id)])
  write_csv_safe(dat, file.path(output_dir, "E2_sample_level_patient_metrics.csv"))
  predictors <- c("distance_score", "model_disagreement", "agreement_score", "sample_level_reliability",
                  "legacy_saved_reliability", "legacy_saved_disagreement")
  set.seed(20261016L); B <- 2000L
  bs <- bind_rows(lapply(predictors, function(v) {
    rho_boot <- replicate(B, { ii <- sample.int(nrow(dat), nrow(dat), TRUE)
      rank_cor(dat[[v]][ii], dat$PCA_ElasticNet_sample_RMSE[ii]) })
    rho_point <- rank_cor(dat[[v]], dat$PCA_ElasticNet_sample_RMSE)
    tibble(score = v, n = sum(complete.cases(dat[[v]], dat$PCA_ElasticNet_sample_RMSE)),
      rho_vs_actual_RMSE = rho_point,
      CI_low = quantile(rho_boot, .025, na.rm = TRUE), CI_high = quantile(rho_boot, .975, na.rm = TRUE),
      bootstrap_B = B, source = if (startsWith(v, "legacy_")) "legacy saved object, provenance not used for primary inference" else "recomputed from canonical split and registered model predictions")
  }))
  write_csv_safe(bs, file.path(output_dir, "E2_sample_level_associations.csv"))
}

run_e2 <- function(x, out, B = 2000L) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  cal_y <- assert_matrix(x$Y_cal_z, "Y_cal_z"); cal_p <- assert_matrix(x$pred_cal_z, "pred_cal_z")
  test_y <- assert_matrix(x$Y_test_z, "Y_test_z"); test_p <- assert_matrix(x$pred_test_z, "pred_test_z")
  train_y <- assert_matrix(x$Y_train_z, "Y_train_z")
  stopifnot(identical(colnames(cal_y), colnames(cal_p)), identical(colnames(test_y), colnames(test_p)))
  cal <- metric_cols(cal_y, cal_p); test <- metric_cols(test_y, test_p)
  native_train <- assert_matrix(x$Y_train_native %||% x$Y_train_z, "Y_train_native (or Y_train_z fallback)")
  scores <- score_variants(cal, train_y, native_train)
  scdf <- tibble(miRNA = cal$miRNA) %>%
    left_join(cal, by = "miRNA", suffix = c("", "_cal")) %>%
    left_join(test, by = "miRNA", suffix = c("_cal", "_test"))
  for (nm in names(scores)) scdf[[nm]] <- unname(scores[[nm]][scdf$miRNA])
  scdf <- scdf %>% mutate(R2_NA = !is.finite(R2_cal), Pearson_NA = !is.finite(Pearson_cal),
                           RMSE_NA = !is.finite(RMSE_cal))
  write_csv_safe(scdf, file.path(out, "E2_miRNA_score_ablation.csv"))
  assoc <- bind_rows(lapply(names(scores), function(s) tibble(
    score = s, n_miRNA = sum(complete.cases(scores[[s]], test$RMSE)),
    rho_test_RMSE = rank_cor(scores[[s]], test$RMSE),
    rho_test_R2 = rank_cor(scores[[s]], test$R2),
    rho_test_Pearson = rank_cor(scores[[s]], test$Pearson))))
  write_csv_safe(assoc, file.path(out, "E2_score_error_associations_point.csv"))
  boot_assoc <- bootstrap_assoc(cal_y, cal_p, test_y, test_p, native_train, B)
  write_csv_safe(boot_assoc, file.path(out, "E2_score_error_associations_bootstrap_raw.csv"))
  ba <- bind_rows(lapply(setdiff(names(boot_assoc), "replicate"), function(v) tibble(
    association = v, estimate_median = median(boot_assoc[[v]], na.rm = TRUE),
    CI_low = quantile(boot_assoc[[v]], .025, na.rm = TRUE),
    CI_high = quantile(boot_assoc[[v]], .975, na.rm = TRUE),
    B_valid = sum(is.finite(boot_assoc[[v]]), na.rm = TRUE), B_requested = B)))
  write_csv_safe(ba, file.path(out, "E2_score_error_associations_bootstrap_summary.csv"))
  cutoffs <- seq(.1, .9, by = .05)
  th <- bind_rows(lapply(names(scores), function(s) {
    if (s != "Rel_equal") return(tibble())
    threshold_metrics(scores[[s]], test, cutoffs, function(r2, p) r2 > 0 & p >= .5) %>% mutate(score = s)
  }))
  write_csv_safe(th, file.path(out, "E2_Rel_threshold_coverage_FNR_point.csv"))
  sens <- bind_rows(lapply(c(.3, .5, .7), function(rmin) {
    threshold_metrics(scores$Rel_equal, test, c(.5, .6, .7),
      function(r2, p) r2 > 0 & p >= rmin) %>% mutate(Pearson_label_threshold = rmin)
  }))
  write_csv_safe(sens, file.path(out, "E2_FNR_label_sensitivity_point.csv"))
  bootfnr <- bootstrap_fnr(cal_y, cal_p, test_y, test_p, B, seq(.1, .9, by = .05), c(.3, .5, .7))
  write_csv_safe(bootfnr, file.path(out, "E2_FNR_bootstrap_raw.csv"))
  write_csv_safe(summarise_boot(bootfnr, c("threshold", "Pearson_label_threshold")),
                 file.path(out, "E2_FNR_bootstrap_summary.csv"))
  # Equal-size score ranking comparisons at 33/74/177 and 3.9/8.7/20.9%.
  sizes <- unique(pmin(c(33L, 74L, 177L), sum(is.finite(scores$Rel_equal))))
  cov <- c(.039, .087, .209)
  for (n in pmax(1L, round(cov * sum(is.finite(scores$Rel_equal))))) sizes <- c(sizes, n)
  sizes <- sort(unique(sizes))
  comparison <- bind_rows(lapply(sizes, function(k) bind_rows(lapply(names(scores), function(s) {
    ord <- rank_score(scores[[s]], names(scores[[s]])); chosen <- head(ord, k)
    mm <- test[match(chosen, test$miRNA), ]
    tibble(target_n = k, score = s, n_selected = length(chosen), coverage = length(chosen)/nrow(test),
      median_R2 = median(mm$R2, na.rm = TRUE), median_Pearson = median(mm$Pearson, na.rm = TRUE),
      median_RMSE = median(mm$RMSE, na.rm = TRUE), rho_test_RMSE = rank_cor(scores[[s]], test$RMSE))
  }))))
  write_csv_safe(comparison, file.path(out, "E2_equal_coverage_score_comparison.csv"))
  sample_level_analysis(x, out)
  oof_path <- "new/results/E1/E1_oof_error_baseline.csv"
  write_csv_safe(tibble(status = "pending_split_matched_source",
    expected_source = oof_path,
    note = "Await split-matched OOF scores; never join pooled E1 cross-split OOF with the fixed E0 test set."),
    file.path(out, "E2_OOF_baseline_status.csv"))
}

bootstrap_external <- function(real, pred, score, B = 2000L, seed = 20261017L) {
  set.seed(seed); n <- nrow(real); out <- vector("list", B)
  mi <- colnames(real)
  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    m <- metric_cols(real[ii, , drop = FALSE], pred[ii, , drop = FALSE])
    ix <- match(mi, names(score))
    out[[b]] <- tibble(replicate = b,
      rho_internal_score_external_RMSE = rank_cor(score[ix], m$RMSE),
      rho_internal_score_external_R2 = rank_cor(score[ix], m$R2),
      rho_internal_score_external_Pearson = rank_cor(score[ix], m$Pearson),
      n_score_RMSE = sum(complete.cases(score[ix], m$RMSE)),
      n_score_R2 = sum(complete.cases(score[ix], m$R2)),
      n_score_Pearson = sum(complete.cases(score[ix], m$Pearson)),
      median_RMSE = median(m$RMSE, na.rm = TRUE), median_R2 = median(m$R2, na.rm = TRUE),
      median_Pearson = median(m$Pearson, na.rm = TRUE))
  }
  bind_rows(out)
}

run_e4 <- function(x, out, B = 2000L, external_rda = "data/livergabage_cancer_data/chol/data.rda") {
  cat("E4 input:", external_rda, "; output:", out, "; bootstrap B:", B, "\n")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  obj <- x$metadata_lightweight_object
  if (is.null(obj) && !is.null(x$lightweight_object_path) && file.exists(x$lightweight_object_path))
    obj <- readRDS(x$lightweight_object_path)
  model_path <- x$metadata_model_path
  if (is.null(obj) || is.null(model_path) || !file.exists(model_path)) {
    write_csv_safe(tibble(status = "pending", reason = "canonical metadata_lightweight_object/model_path absent"),
                   file.path(out, "E4_status.csv")); return(invisible(NULL))
  }
  # Reuse existing, read-only external preprocessing/prediction helpers.
  source("code/run_external_generalization.R", local = TRUE)
  ext <- load_external_data_rda(external_rda)
  built <- build_external_matrices(ext$mrna, ext$mirna, obj)
  model <- readRDS(model_path)
  pred_all <- predict_external_with_model(model, obj, built$x_private_z, built$x_private_all_z)
  common <- intersect(colnames(pred_all), colnames(built$y_true_z))
  real <- built$y_true_z[, common, drop = FALSE]
  pred <- pred_all[, common, drop = FALSE]
  metrics <- metric_cols(real, pred)
  # Fixed internal score only; no external outcome selects weights or features.
  st <- x$score_table
  if (is.null(st) || !all(c("miRNA", "Rel") %in% names(st))) stop("canonical score_table must contain miRNA and Rel")
  score <- setNames(st$Rel, st$miRNA)
  metrics$internal_Rel <- unname(score[metrics$miRNA])
  metrics$internal_R2 <- unname(setNames(st$R2, st$miRNA)[metrics$miRNA])
  metrics$internal_Pearson <- unname(setNames(st$Pearson, st$miRNA)[metrics$miRNA])
  metrics$internal_RMSE <- unname(setNames(st$RMSE, st$miRNA)[metrics$miRNA])
  write_csv_safe(metrics, file.path(out, "E4_external_metrics_by_miRNA.csv"))
  write_csv_safe(tibble(n_overlap_miRNA = nrow(metrics), R2_valid_n = sum(is.finite(metrics$R2)),
    Pearson_valid_n = sum(is.finite(metrics$Pearson)), RMSE_valid_n = sum(is.finite(metrics$RMSE)),
    mean_R2 = mean(metrics$R2, na.rm = TRUE), median_R2 = median(metrics$R2, na.rm = TRUE),
    mean_Pearson = mean(metrics$Pearson, na.rm = TRUE), median_Pearson = median(metrics$Pearson, na.rm = TRUE),
    mean_RMSE = mean(metrics$RMSE, na.rm = TRUE), median_RMSE = median(metrics$RMSE, na.rm = TRUE),
    reconstructable_n = sum(metrics$R2 > 0 & metrics$Pearson >= .5, na.rm = TRUE)),
    file.path(out, "E4_external_overall_metrics.csv"))
  # Coverage and performance of prespecified internal score thresholds on the
  # 742 overlapping miRNAs. External outcomes are evaluation-only.
  rel_sets <- lapply(c(.5, .6, .7), function(cut) st$miRNA[is.finite(st$Rel) & st$Rel >= cut])
  names(rel_sets) <- c("Rel_ge_0.5", "Rel_ge_0.6", "Rel_ge_0.7_canonical_33")
  legacy_sets <- obj$miRNA_sets %||% list()
  if (!is.null(legacy_sets$Rel_ge_0_7)) rel_sets$legacy_saved_33 <- legacy_sets$Rel_ge_0_7
  threshold_summary <- bind_rows(lapply(names(rel_sets), function(nm) {
    selected <- intersect(rel_sets[[nm]], metrics$miRNA)
    mm <- metrics[match(selected, metrics$miRNA), ]
    tibble(selection = nm, internal_selected = length(rel_sets[[nm]]),
      external_overlap = length(selected), overlap_fraction = length(selected)/length(rel_sets[[nm]]),
      metric_valid_n = sum(is.finite(mm$R2) & is.finite(mm$Pearson)),
      R2_valid_n = sum(is.finite(mm$R2)), Pearson_valid_n = sum(is.finite(mm$Pearson)),
      median_R2 = median(mm$R2, na.rm = TRUE), median_Pearson = median(mm$Pearson, na.rm = TRUE),
      median_RMSE = median(mm$RMSE, na.rm = TRUE),
      reconstructable_n = sum(mm$R2 > 0 & mm$Pearson >= .5, na.rm = TRUE),
      reconstructable_fraction = mean(mm$R2 > 0 & mm$Pearson >= .5, na.rm = TRUE))
  }))
  write_csv_safe(threshold_summary, file.path(out, "E4_Rel_threshold_external_coverage.csv"))
  write_csv_safe(bind_rows(lapply(names(rel_sets), function(nm) tibble(selection = nm,
    miRNA = rel_sets[[nm]], shared_with_external = rel_sets[[nm]] %in% metrics$miRNA))),
    file.path(out, "E4_Rel_threshold_feature_lists.csv"))
  # Equal-size selection comparisons. Rankings are frozen from internal data;
  # random sets provide size-matched reference distributions.
  common_names <- metrics$miRNA
  train_native <- as.matrix(x$Y_train_native)
  train_mean <- colMeans(train_native, na.rm = TRUE)
  rank_vectors <- list(
    Internal_Rel = setNames(st$Rel, st$miRNA),
    Internal_Pearson = setNames(st$Pearson, st$miRNA),
    Internal_RMSE = setNames(-st$RMSE, st$miRNA),
    Training_native_mean = setNames(train_mean, colnames(train_native)))
  target_sizes <- sort(unique(pmin(c(33L, 74L, 177L), length(common_names))))
  set.seed(20261018L)
  point_rows <- list(); feature_rows <- list()
  for (k in target_sizes) {
    for (method in names(rank_vectors)) {
      ord <- rank_score(rank_vectors[[method]], intersect(common_names, names(rank_vectors[[method]])))
      chosen <- head(ord, k)
      mm <- metrics[match(chosen, metrics$miRNA), ]
      point_rows[[length(point_rows)+1L]] <- tibble(target_n = k, method = method, replicate = 0L,
        n_selected = length(chosen), median_R2 = median(mm$R2, na.rm = TRUE),
        median_Pearson = median(mm$Pearson, na.rm = TRUE), median_RMSE = median(mm$RMSE, na.rm = TRUE),
        reconstructable_fraction = mean(mm$R2 > 0 & mm$Pearson >= .5, na.rm = TRUE))
      feature_rows[[length(feature_rows)+1L]] <- tibble(target_n = k, method = method, miRNA = chosen)
    }
    for (b in seq_len(200L)) {
      chosen <- sample(common_names, k, replace = FALSE)
      mm <- metrics[match(chosen, metrics$miRNA), ]
      point_rows[[length(point_rows)+1L]] <- tibble(target_n = k, method = "Random", replicate = b,
        n_selected = length(chosen), median_R2 = median(mm$R2, na.rm = TRUE),
        median_Pearson = median(mm$Pearson, na.rm = TRUE), median_RMSE = median(mm$RMSE, na.rm = TRUE),
        reconstructable_fraction = mean(mm$R2 > 0 & mm$Pearson >= .5, na.rm = TRUE))
      feature_rows[[length(feature_rows)+1L]] <- tibble(target_n = k, method = paste0("Random_", b), miRNA = chosen)
    }
  }
  point_df <- bind_rows(point_rows); features_df <- bind_rows(feature_rows)
  write_csv_safe(point_df, file.path(out, "E4_equal_size_selection_raw.csv"))
  write_csv_safe(features_df, file.path(out, "E4_equal_size_selection_feature_lists.csv"))
  write_csv_safe(point_df %>% group_by(target_n, method) %>% summarise(
    n_sets = n(), median_of_set_median_R2 = median(median_R2, na.rm = TRUE),
    median_of_set_median_Pearson = median(median_Pearson, na.rm = TRUE),
    median_of_set_median_RMSE = median(median_RMSE, na.rm = TRUE),
    subset_q025_R2 = quantile(median_R2, .025, na.rm = TRUE),
    subset_q975_R2 = quantile(median_R2, .975, na.rm = TRUE),
    subset_q025_Pearson = quantile(median_Pearson, .025, na.rm = TRUE),
    subset_q975_Pearson = quantile(median_Pearson, .975, na.rm = TRUE),
    subset_q025_RMSE = quantile(median_RMSE, .025, na.rm = TRUE),
    subset_q975_RMSE = quantile(median_RMSE, .975, na.rm = TRUE),
    reconstructable_fraction_median = median(reconstructable_fraction, na.rm = TRUE), .groups = "drop"),
    file.path(out, "E4_equal_size_selection_summary.csv"))
  historical <- readr::read_csv("output/external_generalization/external_metrics_by_miRNA.csv", show_col_types = FALSE) %>%
    filter(model == "PCA_ElasticNet")
  if (nrow(historical)) {
    chk <- metrics %>% select(miRNA, R2, Pearson, RMSE) %>% inner_join(
      historical %>% select(miRNA, old_R2 = R2, old_Pearson = Pearson, old_RMSE = RMSE), by = "miRNA") %>%
      mutate(delta_R2 = R2 - old_R2, delta_Pearson = Pearson - old_Pearson, delta_RMSE = RMSE - old_RMSE)
    write_csv_safe(chk, file.path(out, "E4_historical_prediction_reproduction_check.csv"))
  }
  bs <- bootstrap_external(real, pred, score, B)
  write_csv_safe(bs, file.path(out, "E4_external_patient_bootstrap_raw.csv"))
  boot_summary <- bind_rows(lapply(setdiff(names(bs), "replicate"), function(nm) tibble(
    metric = nm, estimate_median = median(bs[[nm]], na.rm = TRUE),
    CI_low = quantile(bs[[nm]], .025, na.rm = TRUE), CI_high = quantile(bs[[nm]], .975, na.rm = TRUE),
    B_valid = sum(is.finite(bs[[nm]])), B_requested = B,
    effective_miRNA_n_median = if (nm == "rho_internal_score_external_RMSE") median(bs$n_score_RMSE, na.rm = TRUE) else
      if (nm == "rho_internal_score_external_R2") median(bs$n_score_R2, na.rm = TRUE) else
      if (nm == "rho_internal_score_external_Pearson") median(bs$n_score_Pearson, na.rm = TRUE) else NA_real_,
    effective_miRNA_n_min = if (nm == "rho_internal_score_external_RMSE") min(bs$n_score_RMSE, na.rm = TRUE) else
      if (nm == "rho_internal_score_external_R2") min(bs$n_score_R2, na.rm = TRUE) else
      if (nm == "rho_internal_score_external_Pearson") min(bs$n_score_Pearson, na.rm = TRUE) else NA_real_,
    effective_miRNA_n_max = if (nm == "rho_internal_score_external_RMSE") max(bs$n_score_RMSE, na.rm = TRUE) else
      if (nm == "rho_internal_score_external_R2") max(bs$n_score_R2, na.rm = TRUE) else
      if (nm == "rho_internal_score_external_Pearson") max(bs$n_score_Pearson, na.rm = TRUE) else NA_real_)))
  write_csv_safe(boot_summary, file.path(out, "E4_external_patient_bootstrap_summary.csv"))
  # Same-model internal coverage mask: set features absent from external input to
  # the training-centered value (z=0), then quantify the prediction change.
  missing <- setdiff(obj$top_genes, colnames(as.data.frame(ext$mrna)[, setdiff(names(ext$mrna), "patient_id"), drop = FALSE]))
  if (length(missing) && all(obj$top_genes %in% colnames(x$X_test_z))) {
    xx <- x$X_test_z[, obj$top_genes, drop = FALSE]
    xx_mask <- xx
    missing_in_model <- intersect(missing, colnames(xx_mask))
    fill_raw <- obj$x_median[missing_in_model]
    fill_raw[!is.finite(fill_raw)] <- 0
    fill_z <- (fill_raw - obj$x_center[missing_in_model]) / obj$x_scale[missing_in_model]
    fill_z[!is.finite(fill_z)] <- 0
    xx_mask[, missing_in_model] <- matrix(fill_z, nrow(xx_mask), length(fill_z), byrow = TRUE)
    model_internal <- readRDS(model_path)
    pc <- predict(model_internal$pca, newdata = xx[, model_internal$pca$rotation |> rownames(), drop = FALSE])
    pm <- predict(model_internal$pca, newdata = xx_mask[, model_internal$pca$rotation |> rownames(), drop = FALSE])
    kk <- seq_len(model_internal$best_k)
    pf <- predict(model_internal$fit, newx = pc[, kk, drop = FALSE], s = "lambda.min")
    pmm <- predict(model_internal$fit, newx = pm[, kk, drop = FALSE], s = "lambda.min")
    pf <- to_prediction_matrix(pf, rownames(x$Y_test_z), colnames(x$Y_test_z))
    pmm <- to_prediction_matrix(pmm, rownames(x$Y_test_z), colnames(x$Y_test_z))
    fullm <- metric_cols(x$Y_test_z, pf); maskm <- metric_cols(x$Y_test_z, pmm)
    write_csv_safe(tibble(miRNA = fullm$miRNA, full_R2 = fullm$R2, masked_R2 = maskm$R2,
      full_Pearson = fullm$Pearson, masked_Pearson = maskm$Pearson,
      full_RMSE = fullm$RMSE, masked_RMSE = maskm$RMSE,
      missing_gene_count = length(missing), missing_gene_fraction = length(missing)/length(obj$top_genes)),
      file.path(out, "E4_internal_missing_gene_mask_diagnostic.csv"))
    write_csv_safe(tibble(missing_gene = missing), file.path(out, "E4_external_missing_model_genes.csv"))
  } else {
    write_csv_safe(tibble(status = if (length(missing) == 0L) "not_applicable_all_model_genes_present" else "pending", reason = if (length(missing) == 0L) "no missing model genes to mask" else "gene namespace mismatch or canonical X_test_z not aligned with model top_genes"),
                   file.path(out, "E4_internal_missing_gene_mask_status.csv"))
  }
  write_csv_safe(tibble(n_patients = nrow(real), n_overlap_miRNA = ncol(real),
    external_gene_coverage = built$gene_coverage, selected_model = "PCA_ElasticNet",
    model_source = model_path, external_input = external_rda,
    clinical_survival = "not assessed; expression-only external evaluation"),
    file.path(out, "E4_run_metadata.csv"))
  cat("E4 completed:", nrow(real), "patients;", ncol(real), "shared miRNAs; input gene coverage", built$gene_coverage, "\n")
}

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1L]] }
canonical_path <- parse_arg("--canonical", "new/results/E0/canonical_input.rds")
section <- parse_arg("--section", "all")
B <- as.integer(parse_arg("--bootstrap", "2000"))
out_root <- parse_arg("--output-root", "new/results")
external_rda <- parse_arg("--external-input", "data/livergabage_cancer_data/chol/data.rda")
if (!file.exists(canonical_path)) stop("canonical E0 input not found: ", canonical_path)
x <- readRDS(canonical_path)
required <- c("X_train_z", "X_cal_z", "X_test_z", "Y_train_z", "Y_cal_z", "Y_test_z",
              "pred_cal_z", "pred_test_z", "miRNA_names", "gene_names", "score_table")
missing <- setdiff(required, names(x))
if (length(missing)) stop("canonical input missing fields: ", paste(missing, collapse = ", "))
if (section %in% c("E2", "all")) run_e2(x, file.path(out_root, "E2"), B)
if (section == "sample") sample_level_analysis(x, file.path(out_root, "E2"))
if (section %in% c("E4", "all")) run_e4(x, file.path(out_root, "E4"), B, external_rda)
