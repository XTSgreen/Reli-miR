#!/usr/bin/env Rscript
# Completes paired uncertainty estimates omitted from the first E2/E4 run.
# All bootstrap intervals are conditional on fitted models and use paired
# patient resampling. Scores are re-estimated within each E2 draw.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

fast_metrics <- function(real, pred) {
  real <- as.matrix(real); pred <- as.matrix(pred)
  stopifnot(identical(dim(real), dim(pred)))
  if (is.null(colnames(real))) colnames(real) <- paste0("v", seq_len(ncol(real)))
  if (is.null(colnames(pred))) colnames(pred) <- colnames(real)
  stopifnot(identical(colnames(real), colnames(pred)))
  y_valid <- is.finite(real)
  pair_valid <- y_valid & is.finite(pred)
  y0 <- real; y0[!y_valid] <- 0
  yp <- real; pp <- pred
  yp[!pair_valid] <- 0; pp[!pair_valid] <- 0
  n_y <- colSums(y_valid); n_pair <- colSums(pair_valid)
  sy <- colSums(y0); sy2 <- colSums(y0 * y0)
  sp <- colSums(pp); sp2 <- colSums(pp * pp); syp <- colSums(yp * pp)
  # Compute R2 from the same pairwise-complete vector used by metric_cols;
  # direct centering also avoids cancellation for near-constant resamples.
  sst <- numeric(ncol(real))
  for (j in seq_len(ncol(real))) {
    ok <- pair_valid[, j]
    if (sum(ok) >= 3L) sst[j] <- sum((real[ok, j] - mean(real[ok, j]))^2)
  }
  vy_pair <- sy2 # replaced below with pair-only sufficient statistics
  yp2 <- colSums(yp * yp)
  vp <- sp2 - ifelse(n_pair > 0, sp^2 / n_pair, 0)
  vy_pair <- yp2 - ifelse(n_pair > 0, colSums(yp)^2 / n_pair, 0)
  cov_pair <- syp - ifelse(n_pair > 0, colSums(yp) * sp / n_pair, 0)
  rss <- colSums((yp - pp)^2)
  r2 <- rep(NA_real_, ncol(real)); ok_r2 <- n_pair >= 3 & is.finite(sst) & sst > 0
  r2[ok_r2] <- 1 - rss[ok_r2] / sst[ok_r2]
  pearson <- rep(NA_real_, ncol(real))
  # Match metric_cols: the nonconstant check uses each complete vector, then
  # correlation itself uses pairwise-complete observations.
  sd_y_all <- apply(real, 2, sd, na.rm = TRUE)
  sd_p_all <- apply(pred, 2, sd, na.rm = TRUE)
  ok_cor <- n_pair >= 2 & is.finite(sd_y_all) & is.finite(sd_p_all) & sd_y_all > 1e-8 & sd_p_all > 1e-8 & vy_pair > 0 & vp > 0
  pearson[ok_cor] <- cov_pair[ok_cor] / sqrt(vy_pair[ok_cor] * vp[ok_cor])
  rmse <- rep(NA_real_, ncol(real)); ok_rmse <- n_pair > 0
  rmse[ok_rmse] <- sqrt(rss[ok_rmse] / n_pair[ok_rmse])
  tibble(miRNA = colnames(real), R2 = r2, Pearson = pearson, RMSE = rmse,
         n_real = n_y, n_paired = n_pair)
}

slow_metrics <- function(real, pred) {
  real <- as.matrix(real); pred <- as.matrix(pred)
  r2 <- pearson <- rmse <- rep(NA_real_, ncol(real))
  for (j in seq_len(ncol(real))) {
    ok <- is.finite(real[, j]) & is.finite(pred[, j]); a <- real[ok, j]; b <- pred[ok, j]
    if (length(a) >= 3L) {
      tss <- sum((a - mean(a))^2)
      if (tss > 0) r2[j] <- 1 - sum((a - b)^2) / tss
      if (sd(a) > 0 && sd(b) > 0) pearson[j] <- suppressWarnings(cor(a, b))
      rmse[j] <- sqrt(mean((a - b)^2))
    }
  }
  tibble(miRNA = colnames(real), R2 = r2, Pearson = pearson, RMSE = rmse)
}

rho <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}

score_set <- function(cal_m, train_native) {
  r2 <- setNames(cal_m$R2, cal_m$miRNA)
  pr <- setNames(cal_m$Pearson, cal_m$miRNA)
  rm <- setNames(cal_m$RMSE, cal_m$miRNA)
  native_mean <- colMeans(train_native, na.rm = TRUE)
  native_var <- apply(train_native, 2, var, na.rm = TRUE)
  r2c <- pmin(pmax(r2, 0), 1); prc <- pmin(pmax(pr, 0), 1)
  list(
    Rel_equal = .5 * r2c + .5 * prc,
    R2_raw = r2,
    R2_clipped = r2c,
    Pearson = pr,
    Pearson_clipped = prc,
    Rel_R2_weight_025 = .25 * r2c + .75 * prc,
    Rel_R2_weight_075 = .75 * r2c + .25 * prc,
    Calibration_RMSE_rank = -rm,
    Training_native_mean = native_mean,
    Training_native_variance = native_var
  )
}

summarise_ci <- function(dat, keys, value) {
  dat %>% group_by(across(all_of(keys))) %>% summarise(
    B_valid = sum(is.finite(.data[[value]])),
    bootstrap_median = median(.data[[value]], na.rm = TRUE),
    CI_q025 = quantile(.data[[value]], .025, na.rm = TRUE, names = FALSE),
    CI_q975 = quantile(.data[[value]], .975, na.rm = TRUE, names = FALSE),
    .groups = "drop")
}

paired_metrics <- function(real, pred_list, idx) {
  bind_rows(lapply(names(pred_list), function(nm) {
    m <- fast_metrics(real[idx, , drop = FALSE], pred_list[[nm]][idx, , drop = FALSE])
    tibble(model = nm, mean_R2 = mean(m$R2, na.rm = TRUE),
      median_Pearson = median(m$Pearson, na.rm = TRUE),
      mean_RMSE = mean(m$RMSE, na.rm = TRUE),
      n_R2 = sum(is.finite(m$R2)), n_Pearson = sum(is.finite(m$Pearson)),
      n_RMSE = sum(is.finite(m$RMSE)))
  }))
}

run_model_bootstrap <- function(x, B, out_dir) {
  y <- as.matrix(x$Y_test_z)
  models <- c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")
  pred <- setNames(lapply(models, function(m) readRDS(file.path("output/predictions", paste0(m, "_prediction_test_z.rds")))), models)
  for (nm in models) {
    pred[[nm]] <- as.matrix(pred[[nm]])
    if (!identical(dim(pred[[nm]]), dim(y)) || !identical(dimnames(pred[[nm]]), dimnames(y)))
      stop("prediction misalignment for ", nm)
  }
  point <- paired_metrics(y, pred, seq_len(nrow(y))) %>% mutate(bootstrap_scope = "fixed fitted model; conditional on canonical E0 split")
  set.seed(20261019L)
  raw <- vector("list", B)
  for (b in seq_len(B)) {
    ii <- sample.int(nrow(y), nrow(y), replace = TRUE)
    raw[[b]] <- paired_metrics(y, pred, ii) %>% mutate(replicate = b)
  }
  boot <- bind_rows(raw)
  write_csv(point, file.path(out_dir, "E2_complete_model_metrics_point.csv"))
  write_csv(boot, file.path(out_dir, "E2_complete_model_bootstrap_raw.csv"))
  metrics <- c("mean_R2", "median_Pearson", "mean_RMSE")
  by_model <- bind_rows(lapply(metrics, function(m) {
    s <- summarise_ci(boot, "model", m)
    p <- point %>% select(model, point_estimate = all_of(m))
    left_join(s, p, by = "model") %>% mutate(metric = m,
      bootstrap_scope = "paired patient bootstrap; fixed fitted models; conditional on canonical E0 split")
  }))
  write_csv(by_model, file.path(out_dir, "E2_complete_model_bootstrap_summary.csv"))
  pairs <- combn(models, 2, simplify = FALSE)
  contrasts <- bind_rows(lapply(pairs, function(pair) bind_rows(lapply(metrics, function(m) {
    a <- boot %>% filter(model == pair[1]) %>% arrange(replicate) %>% pull(all_of(m))
    b <- boot %>% filter(model == pair[2]) %>% arrange(replicate) %>% pull(all_of(m))
    p <- point %>% filter(model %in% pair) %>% arrange(match(model, pair)) %>% pull(all_of(m))
    d <- a - b
    tibble(model_a = pair[1], model_b = pair[2], metric = m,
      contrast_definition = paste0(pair[1], " minus ", pair[2]),
      point_difference = p[1] - p[2], bootstrap_median_difference = median(d, na.rm = TRUE),
      CI_q025 = quantile(d, .025, na.rm = TRUE), CI_q975 = quantile(d, .975, na.rm = TRUE),
      B_valid = sum(is.finite(d)), B_requested = B,
      bootstrap_scope = "paired patient bootstrap; fixed fitted models; conditional on canonical E0 split")
  }))))
  write_csv(contrasts, file.path(out_dir, "E2_complete_model_paired_contrasts.csv"))
}

run_score_bootstrap <- function(x, B, out_dir) {
  ycal <- as.matrix(x$Y_cal_z); pcal <- as.matrix(x$pred_cal_z)
  ytest <- as.matrix(x$Y_test_z); ptest <- as.matrix(x$pred_test_z)
  ytrain <- as.matrix(x$Y_train_native)
  vars <- score_set(fast_metrics(ycal, pcal), ytrain)
  testpoint <- fast_metrics(ytest, ptest)
  point_rho <- vapply(vars, function(s) rho(s, testpoint$RMSE), numeric(1))
  rel_point <- unname(point_rho[["Rel_equal"]])
  point <- tibble(score = names(point_rho), point_rho_vs_test_RMSE = as.numeric(point_rho),
    point_rho_baseline_minus_Rel = as.numeric(point_rho - rel_point),
    contrast_direction = "baseline rho minus Rel rho; positive supports stronger negative Rel association",
    score_estimation = "calibration and training summaries on fixed canonical samples")
  write_csv(point, file.path(out_dir, "E2_complete_score_contrasts_point.csv"))
  # Confirm score-point agreement with the E2 point table before resampling.
  existing <- read_csv("new/results/E2/E2_score_error_associations_point.csv", show_col_types = FALSE)
  check <- point %>% inner_join(existing %>% select(score, existing_rho = rho_test_RMSE), by = "score") %>%
    mutate(abs_point_delta = abs(point_rho_vs_test_RMSE - existing_rho))
  write_csv(check, file.path(out_dir, "E2_complete_score_point_reproduction_check.csv"))
  if (nrow(check) < 7L || any(check$abs_point_delta > 1e-8, na.rm = TRUE))
    stop("E2 point score associations do not reproduce the existing E2 table")
  set.seed(20261020L)
  raw <- vector("list", B)
  for (b in seq_len(B)) {
    itrain <- sample.int(nrow(ytrain), nrow(ytrain), replace = TRUE)
    ical <- sample.int(nrow(ycal), nrow(ycal), replace = TRUE)
    itest <- sample.int(nrow(ytest), nrow(ytest), replace = TRUE)
    train_b <- ytrain[itrain, , drop = FALSE]
    cal_b <- fast_metrics(ycal[ical, , drop = FALSE], pcal[ical, , drop = FALSE])
    test_b <- fast_metrics(ytest[itest, , drop = FALSE], ptest[itest, , drop = FALSE])
    sb <- score_set(cal_b, train_b)
    rh <- vapply(sb, function(s) rho(s, test_b$RMSE), numeric(1))
    rel <- unname(rh[["Rel_equal"]])
    raw[[b]] <- tibble(replicate = b, score = names(rh), rho_vs_test_RMSE = as.numeric(rh),
      rho_baseline_minus_Rel = as.numeric(rh - rel), n_effective_miRNA = vapply(sb, function(s)
        sum(is.finite(s) & is.finite(test_b$RMSE)), integer(1)),
      score_estimation = "re-estimated from jointly resampled train-native, calibration and test patients",
      bootstrap_scope = "paired by replicate across all scoring rules; fixed fitted prediction matrices")
  }
  boot <- bind_rows(raw)
  write_csv(boot, file.path(out_dir, "E2_complete_score_joint_bootstrap_raw.csv"))
  intervals <- bind_rows(lapply(unique(boot$score), function(s) {
    z <- boot %>% filter(score == s)
    tibble(score = s, point_rho_vs_test_RMSE = point_rho[[s]],
      bootstrap_median_rho = median(z$rho_vs_test_RMSE, na.rm = TRUE),
      rho_CI_q025 = quantile(z$rho_vs_test_RMSE, .025, na.rm = TRUE),
      rho_CI_q975 = quantile(z$rho_vs_test_RMSE, .975, na.rm = TRUE),
      point_rho_baseline_minus_Rel = point_rho[[s]] - rel_point,
      bootstrap_median_baseline_minus_Rel = median(z$rho_baseline_minus_Rel, na.rm = TRUE),
      contrast_CI_q025 = quantile(z$rho_baseline_minus_Rel, .025, na.rm = TRUE),
      contrast_CI_q975 = quantile(z$rho_baseline_minus_Rel, .975, na.rm = TRUE),
      B_valid = sum(is.finite(z$rho_baseline_minus_Rel)), B_requested = B,
      effective_miRNA_n_median = median(z$n_effective_miRNA, na.rm = TRUE),
      score_estimation = first(z$score_estimation), bootstrap_scope = first(z$bootstrap_scope))
  }))
  write_csv(intervals, file.path(out_dir, "E2_complete_score_joint_bootstrap_summary.csv"))
}

run_sample_bootstrap <- function(B, out_dir) {
  dat <- read_csv(file.path(out_dir, "E2_sample_level_patient_metrics.csv"), show_col_types = FALSE)
  outcome <- dat$PCA_ElasticNet_sample_RMSE
  predictors <- c("distance_score", "model_disagreement", "agreement_score", "sample_level_reliability",
                  "legacy_saved_reliability", "legacy_saved_disagreement")
  set.seed(20261016L)
  raw <- bind_rows(lapply(predictors, function(v) {
    z <- dat[[v]]
    rboot <- nboot <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      ii <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
      rboot[b] <- rho(z[ii], outcome[ii])
      nboot[b] <- sum(complete.cases(z[ii], outcome[ii]))
    }
    rho_point <- rho(z, outcome)
    tibble(score = v, replicate = seq_len(B), rho_vs_actual_RMSE = as.numeric(rboot),
      rho_point = rho_point, n_effective = sum(complete.cases(z, outcome)),
      n_boot_effective = nboot,
      bootstrap_scope = "paired resampling of 72 test patients; fixed sample-level score and canonical PCA_ElasticNet error")
  }))
  write_csv(raw, file.path(out_dir, "E2_complete_sample_bootstrap_raw.csv"))
  summary <- raw %>% group_by(score) %>% summarise(
    point_rho_vs_actual_RMSE = first(rho_point),
    bootstrap_median_rho = median(rho_vs_actual_RMSE, na.rm = TRUE),
    CI_q025 = quantile(rho_vs_actual_RMSE, .025, na.rm = TRUE),
    CI_q975 = quantile(rho_vs_actual_RMSE, .975, na.rm = TRUE),
    n_effective = first(n_effective), B_valid = sum(is.finite(rho_vs_actual_RMSE)),
    B_requested = B, bootstrap_scope = first(bootstrap_scope), .groups = "drop")
  write_csv(summary, file.path(out_dir, "E2_complete_sample_bootstrap_summary.csv"))
}

run_external_bootstrap <- function(x, B, out_dir, external_rda = "data/livergabage_cancer_data/chol/data.rda") {
  source("code/run_external_generalization.R", local = TRUE)
  obj_path <- x$lightweight_object_path
  model_path <- x$metadata_model_path
  obj <- readRDS(obj_path); model <- readRDS(model_path)
  ext <- load_external_data_rda(external_rda)
  built <- build_external_matrices(ext$mrna, ext$mirna, obj)
  pred_all <- predict_external_with_model(model, obj, built$x_private_z, built$x_private_all_z)
  common <- intersect(colnames(pred_all), colnames(built$y_true_z))
  real <- built$y_true_z[, common, drop = FALSE]
  pred <- pred_all[, common, drop = FALSE]
  stopifnot(nrow(real) == 35L, ncol(real) == 742L,
            identical(rownames(real), rownames(pred)), identical(colnames(real), colnames(pred)))
  internal <- setNames(x$score_table$Rel, x$score_table$miRNA)[common]
  set.seed(20261017L)
  raw <- vector("list", B)
  for (b in seq_len(B)) {
    ii <- sample.int(nrow(real), nrow(real), replace = TRUE)
    m <- fast_metrics(real[ii, , drop = FALSE], pred[ii, , drop = FALSE])
    raw[[b]] <- tibble(replicate = b,
      mean_R2 = mean(m$R2, na.rm = TRUE), median_R2 = median(m$R2, na.rm = TRUE),
      mean_Pearson = mean(m$Pearson, na.rm = TRUE), median_Pearson = median(m$Pearson, na.rm = TRUE),
      mean_RMSE = mean(m$RMSE, na.rm = TRUE), median_RMSE = median(m$RMSE, na.rm = TRUE),
      rho_internal_Rel_vs_RMSE = rho(internal[m$miRNA], m$RMSE),
      rho_internal_Rel_vs_R2 = rho(internal[m$miRNA], m$R2),
      rho_internal_Rel_vs_Pearson = rho(internal[m$miRNA], m$Pearson),
      n_R2 = sum(is.finite(m$R2)), n_Pearson = sum(is.finite(m$Pearson)), n_RMSE = sum(is.finite(m$RMSE)),
      bootstrap_scope = "paired resampling of 35 external patients; fixed PCA_ElasticNet model and internal Rel scores")
  }
  boot <- bind_rows(raw)
  previous <- read_csv(file.path(out_dir, "E4_external_patient_bootstrap_raw.csv"), show_col_types = FALSE)
  compare <- boot %>% inner_join(previous, by = "replicate", suffix = c("_fast", "_existing")) %>%
    transmute(replicate,
      delta_median_R2 = median_R2_fast - median_R2_existing,
      delta_median_Pearson = median_Pearson_fast - median_Pearson_existing,
      delta_median_RMSE = median_RMSE_fast - median_RMSE_existing,
      delta_rho_RMSE = rho_internal_Rel_vs_RMSE - rho_internal_score_external_RMSE)
  write_csv(compare, file.path(out_dir, "E4_complete_bootstrap_reproduction_check.csv"))
  if (nrow(compare) != B || any(abs(as.matrix(compare[-1])) > 1e-8, na.rm = TRUE))
    stop("E4 fixed-seed external bootstrap does not reproduce existing B=2000 summaries")
  point <- fast_metrics(real, pred)
  point_summary <- tibble(n_patient = nrow(real), n_miRNA = ncol(real),
    point_mean_R2 = mean(point$R2, na.rm = TRUE), point_median_R2 = median(point$R2, na.rm = TRUE),
    point_mean_Pearson = mean(point$Pearson, na.rm = TRUE), point_median_Pearson = median(point$Pearson, na.rm = TRUE),
    point_mean_RMSE = mean(point$RMSE, na.rm = TRUE), point_median_RMSE = median(point$RMSE, na.rm = TRUE),
    point_rho_internal_Rel_vs_RMSE = rho(internal[point$miRNA], point$RMSE),
    point_rho_internal_Rel_vs_R2 = rho(internal[point$miRNA], point$R2),
    point_rho_internal_Rel_vs_Pearson = rho(internal[point$miRNA], point$Pearson),
    model_source = model_path, score_source = "E0 score_table; fixed across external bootstrap", external_input = external_rda)
  write_csv(point_summary, file.path(out_dir, "E4_complete_external_point_summary.csv"))
  write_csv(boot, file.path(out_dir, "E4_complete_external_bootstrap_raw.csv"))
  summaries <- bind_rows(lapply(c("mean_R2", "median_R2", "mean_Pearson", "median_Pearson", "mean_RMSE", "median_RMSE",
                                  "rho_internal_Rel_vs_RMSE", "rho_internal_Rel_vs_R2", "rho_internal_Rel_vs_Pearson"), function(v) {
    tibble(metric = v, point_estimate = point_summary[[paste0("point_", v)]],
      bootstrap_median = median(boot[[v]], na.rm = TRUE),
      CI_q025 = quantile(boot[[v]], .025, na.rm = TRUE), CI_q975 = quantile(boot[[v]], .975, na.rm = TRUE),
      B_valid = sum(is.finite(boot[[v]])), B_requested = B,
      bootstrap_scope = "paired 35-patient bootstrap; fixed PCA_ElasticNet and fixed internal Rel")
  }))
  write_csv(summaries, file.path(out_dir, "E4_complete_external_bootstrap_summary.csv"))
}

sanity_metrics <- function() {
  set.seed(29)
  y <- matrix(rnorm(48), 12, 4, dimnames = list(NULL, paste0("v", 1:4)))
  p <- y + matrix(rnorm(48, sd = .3), 12, 4)
  y[2, 1] <- NA_real_; p[5, 2] <- NA_real_
  y[, 3] <- 2; p[, 3] <- rnorm(12)
  y[, 4] <- rnorm(12); p[, 4] <- 5
  a <- fast_metrics(y, p); b <- slow_metrics(y, p)
  for (metric in c("R2", "Pearson", "RMSE")) {
    av <- a[[metric]]; bv <- b[[metric]]
    if (any(is.finite(av) != is.finite(bv))) {
      message("finite mismatch ", metric, ": fast=", paste(av, collapse = ","),
              "; slow=", paste(bv, collapse = ","),
              "; fast_mask=", paste(is.finite(av), collapse = ","),
              "; slow_mask=", paste(is.finite(bv), collapse = ","))
      stop("finite-mask mismatch in small metric comparison")
    }
    ok <- is.finite(av) & is.finite(bv)
    stopifnot(max(abs(av[ok] - bv[ok])) < 1e-10)
  }
  stopifnot(is.na(a$R2[3]), is.na(a$Pearson[3]), is.na(a$Pearson[4]))
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i >= length(args)) default else args[[i+1L]] }
canonical <- arg("--canonical", "new/results/E0/canonical_input.rds")
B <- as.integer(arg("--bootstrap", "2000"))
section <- arg("--section", "all")
out_root <- arg("--output-root", "new/results")
external_rda <- arg("--external-input", "data/livergabage_cancer_data/chol/data.rda")
if (!file.exists(canonical)) stop("canonical E0 input missing: ", canonical)
sanity_metrics()
x <- readRDS(canonical)
dir.create(file.path(out_root, "E2"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "E4"), recursive = TRUE, showWarnings = FALSE)
if (section %in% c("models", "all")) run_model_bootstrap(x, B, file.path(out_root, "E2"))
if (section %in% c("scores", "all")) run_score_bootstrap(x, B, file.path(out_root, "E2"))
if (section %in% c("sample", "all")) run_sample_bootstrap(B, file.path(out_root, "E2"))
if (section %in% c("external", "all")) run_external_bootstrap(x, B, file.path(out_root, "E4"), external_rda)
cat("sanity_metrics=passed; section=", section, "; bootstrap_B=", B, "\n", sep = "")
