#!/usr/bin/env Rscript
# Fixed-split, paired-permutation negative control for the E3 comparison.
# Each permutation shuffles only the 222 training mRNA/miRNA pairings; calibration
# and test patients and their observed miRNA values remain untouched.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "all"
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
canonical_path <- if (length(args) >= 2L) args[[2L]] else file.path(root, "new", "results", "E0", "canonical_input.rds")
out_root <- if (identical(mode, "pilot")) file.path(root, "new", "results", "_pilot", "E1_shuffled") else file.path(root, "new", "results", "E1_shuffled")
progress_path <- file.path(out_root, "current_progress.csv")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

N_LAMBDA <- 100L                 # glmnet default path length
LAMBDA_MIN_RATIO <- 0.01         # n < p, matching glmnet default
RIDGE_ALPHA_ZERO_EPS <- 0.001   # glmnet alpha=0 lambda_max convention
TOP_MRNA <- 3000L
N_FOLDS <- 5L
PERM_SEEDS <- 20261014:20261033

log_progress <- function(seed, foldstage, status = "running", model = "Shuffled_Ridge", config = NA_integer_, elapsed_seconds = NA_real_) {
  row <- data.frame(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), seed = as.integer(seed),
    model = model, config = as.integer(config), foldstage = as.character(foldstage), status = status, elapsed_seconds = as.numeric(elapsed_seconds),
    check.names = FALSE)
  write.table(row, progress_path, sep = ",", row.names = FALSE, col.names = !file.exists(progress_path),
    append = file.exists(progress_path), na = "NA", qmethod = "double")
}

stratified_folds <- function(group, k = N_FOLDS, seed = 30271014L) {
  set.seed(seed)
  group <- as.character(group); group[is.na(group) | !nzchar(group)] <- "<NA>"
  folds <- vector("list", k)
  for (g in unique(group)) {
    ix <- sample(which(group == g)); assignment <- rep(seq_len(k), length.out = length(ix))
    for (f in seq_len(k)) folds[[f]] <- c(folds[[f]], ix[assignment == f])
  }
  lapply(folds, sort)
}

impute_train <- function(train, other, label) {
  med <- apply(train, 2L, function(z) { z <- z[is.finite(z)]; if (!length(z)) 0 else median(z) })
  for (j in seq_len(ncol(train))) {
    train[!is.finite(train[, j]), j] <- med[j]
    other[!is.finite(other[, j]), j] <- med[j]
  }
  list(train = train, other = other, median = med)
}

fit_x_fold <- function(Xtr_raw, Xva_raw, top_mrna = TOP_MRNA) {
  keep <- colMeans(is.finite(Xtr_raw)) >= 0.8
  Xtr <- Xtr_raw[, keep, drop = FALSE]; Xva <- Xva_raw[, keep, drop = FALSE]
  imp <- impute_train(Xtr, Xva, "X"); Xtr <- imp$train; Xva <- imp$other
  v <- apply(Xtr, 2L, var); keep_v <- is.finite(v) & v > 1e-8
  Xtr <- Xtr[, keep_v, drop = FALSE]; Xva <- Xva[, keep_v, drop = FALSE]; v <- v[keep_v]
  genes <- names(sort(v, decreasing = TRUE))[seq_len(min(top_mrna, length(v)))]
  Xtr <- Xtr[, genes, drop = FALSE]; Xva <- Xva[, genes, drop = FALSE]
  center <- colMeans(Xtr); scale <- apply(Xtr, 2L, sd); scale[!is.finite(scale) | scale < 1e-8] <- 1
  Xtr_z <- sweep(sweep(Xtr, 2L, center, "-"), 2L, scale, "/")
  Xva_z <- sweep(sweep(Xva, 2L, center, "-"), 2L, scale, "/")
  Xtr_z[!is.finite(Xtr_z)] <- 0; Xva_z[!is.finite(Xva_z)] <- 0
  list(train = Xtr_z, valid = Xva_z, genes = genes, median = imp$median[genes], center = center, scale = scale)
}

fit_y_fold <- function(Ytr_raw, Yva_raw) {
  keep <- colMeans(is.finite(Ytr_raw)) >= 0.8
  Ytr <- Ytr_raw[, keep, drop = FALSE]; Yva <- Yva_raw[, keep, drop = FALSE]
  imp <- impute_train(Ytr, Yva, "Y"); Ytr <- imp$train; Yva <- imp$other
  v <- apply(Ytr, 2L, var); keep_v <- is.finite(v) & v > 1e-8
  Ytr <- Ytr[, keep_v, drop = FALSE]; Yva <- Yva[, keep_v, drop = FALSE]
  center <- colMeans(Ytr); scale <- apply(Ytr, 2L, sd); scale[!is.finite(scale) | scale < 1e-8] <- 1
  Ytr_z <- sweep(sweep(Ytr, 2L, center, "-"), 2L, scale, "/")
  Yva_z <- sweep(sweep(Yva, 2L, center, "-"), 2L, scale, "/")
  Ytr_z[!is.finite(Ytr_z)] <- 0; Yva_z[!is.finite(Yva_z)] <- 0
  list(train_z = Ytr_z, valid_z = Yva_z, valid_native = Yva, features = colnames(Ytr),
       median = imp$median[colnames(Ytr)], center = center, scale = scale)
}

lambda_max_mgaussian_ridge <- function(X, Y) {
  # glmnet's alpha=0 path uses alpha=0.001 for lambda_max; the remaining path
  # is geometric. Verified against glmnet on a small multiresponse example.
  grad <- crossprod(X, Y) / nrow(X)
  max(sqrt(rowSums(grad^2))) / RIDGE_ALPHA_ZERO_EPS
}

ridge_spectral_context <- function(Xtr_z, Xva_z, Ytr_z) {
  # Dual/SVD solution for the alpha=0 Frobenius ridge objective used by
  # glmnet(family='mgaussian'). X is already train-centered/scaled, so the
  # intercept is represented by the response mean outside the solve.
  gram <- tcrossprod(Xtr_z)
  eg <- eigen(gram, symmetric = TRUE)
  d2 <- pmax(eg$values, 0)
  U <- eg$vectors
  cross_val_train <- Xva_z %*% t(Xtr_z)
  list(U = U, d2 = d2, XvaU = cross_val_train %*% U,
       q = crossprod(U, Ytr_z), n_train = nrow(Xtr_z))
}

ridge_predict_z <- function(ctx, lambda, path_index) {
  if (path_index == 1L) return(matrix(0, nrow(ctx$XvaU), ncol(ctx$q)))
  # Xnew X' U already contributes the singular-value factor; divide only by
  # (d^2 + n*lambda) here to obtain Xnew V D (D^2+n*lambda)^-1 U'Y.
  shrink <- 1 / (ctx$d2 + ctx$n_train * lambda)
  ctx$XvaU %*% sweep(ctx$q, 1L, shrink, "*")
}

predict_dual_z <- function(Xtr_z, Xnew_z, Ytr_z, lambda, path_index) {
  ctx <- ridge_spectral_context(Xtr_z, Xnew_z, Ytr_z)
  ridge_predict_z(ctx, lambda, path_index)
}

metric_table <- function(real, pred, scale = rep(1, ncol(real))) {
  r <- vapply(seq_len(ncol(real)), function(j) {
    if (sd(real[, j]) < 1e-12 || sd(pred[, j]) < 1e-12) return(NA_real_)
    cor(real[, j], pred[, j])
  }, numeric(1))
  denom <- colSums(sweep(real, 2L, colMeans(real), "-")^2)
  sse <- colSums((real - pred)^2)
  data.frame(miRNA = colnames(real), R2 = 1 - sse / denom, Pearson = r,
    RMSE = sqrt(colMeans((real - pred)^2)),
    RMSE_normalized = sqrt(colMeans((real - pred)^2)) / scale)
}

run_permutation <- function(seed, dat, split, fold_ids, x_fold_cache) {
  target <- file.path(out_root, sprintf("perm_%02d_seed_%d.rds", seed - min(PERM_SEEDS) + 1L, seed))
  summary_target <- file.path(out_root, sprintf("perm_%02d_summary.csv", seed - min(PERM_SEEDS) + 1L))
  if (file.exists(target) && file.exists(summary_target)) {
    old <- tryCatch(readRDS(target), error = function(e) NULL)
    if (!is.null(old) && identical(old$status, "complete") && identical(old$seed, seed)) {
      log_progress(seed, "complete_result_verified_skip", "skipped_complete")
      return(invisible(old))
    }
  }
  t0 <- proc.time()[["elapsed"]]
  train_ids <- split$train; cal_ids <- split$cal; test_ids <- split$test
  Ytr_true <- dat$Y[train_ids, , drop = FALSE]
  set.seed(seed)
  perm <- sample(seq_along(train_ids))
  Ytr_shuffled <- Ytr_true[perm, , drop = FALSE]
  rownames(Ytr_shuffled) <- train_ids
  permutation_map <- data.frame(x_patient_id = train_ids, shuffled_y_source_id = train_ids[perm])
  log_progress(seed, "training_pairing_permuted")

  ratios <- exp(seq(0, log(LAMBDA_MIN_RATIO), length.out = N_LAMBDA))
  fold_cache <- vector("list", length(fold_ids))
  for (fi in seq_along(fold_ids)) {
    va <- fold_ids[[fi]]; fit <- setdiff(seq_along(train_ids), va)
    xfc <- x_fold_cache[[fi]]
    yfc <- fit_y_fold(Ytr_shuffled[fit, , drop = FALSE], Ytr_shuffled[va, , drop = FALSE])
    lmax <- lambda_max_mgaussian_ridge(xfc$train, yfc$train_z)
    lambdas <- lmax * ratios
    spec <- ridge_spectral_context(xfc$train, xfc$valid, yfc$train_z)
    sse <- numeric(N_LAMBDA)
    for (li in seq_len(N_LAMBDA)) {
      pz <- ridge_predict_z(spec, lambdas[li], li)
      sse[li] <- sum((yfc$valid_z - pz)^2)
    }
    fold_cache[[fi]] <- list(valid_indices = va, valid_ids = train_ids[va], yfc = yfc,
      spec = spec, lambdas = lambdas, lambda_max = lmax, sse = sse,
      n_values = length(yfc$valid_z))
    log_progress(seed, paste0("fold_", fi, "_of_", N_FOLDS, "_ridge_path_complete"), elapsed_seconds = proc.time()[["elapsed"]] - t0)
  }
  cv_sse <- Reduce(`+`, lapply(fold_cache, `[[`, "sse"))
  cv_n <- sum(vapply(fold_cache, `[[`, numeric(1), "n_values"))
  best_index <- which.min(cv_sse / cv_n)
  best_ratio <- ratios[best_index]
  oof_pred <- matrix(NA_real_, length(train_ids), ncol(dat$Y), dimnames = list(train_ids, colnames(dat$Y)))
  oof_truth_shuffled <- oof_pred
  oof_norm_sse <- numeric(ncol(dat$Y)); names(oof_norm_sse) <- colnames(dat$Y)
  oof_norm_n <- numeric(ncol(dat$Y)); names(oof_norm_n) <- colnames(dat$Y)
  for (fi in seq_along(fold_cache)) {
    fc <- fold_cache[[fi]]; yfc <- fc$yfc
    pz <- ridge_predict_z(fc$spec, fc$lambdas[best_index], best_index)
    p_native <- sweep(sweep(pz, 2L, yfc$scale, "*"), 2L, yfc$center, "+")
    cols <- yfc$features; rows <- fc$valid_indices
    oof_pred[rows, cols] <- p_native
    oof_truth_shuffled[rows, cols] <- yfc$valid_native
    oof_norm_sse[cols] <- oof_norm_sse[cols] + colSums((yfc$valid_z - pz)^2)
    oof_norm_n[cols] <- oof_norm_n[cols] + nrow(yfc$valid_z)
  }
  oof_rmse_native <- sqrt(colMeans((oof_truth_shuffled - oof_pred)^2, na.rm = TRUE))
  oof_rmse_normalized <- sqrt(oof_norm_sse / oof_norm_n)
  oof_rmse_tbl <- data.frame(miRNA = names(oof_rmse_native),
    OOF_RMSE_native = as.numeric(oof_rmse_native),
    OOF_RMSE_normalized_fold_scale = as.numeric(oof_rmse_normalized))

  # Refit preprocessing and the same relative ridge penalty on all 222 shuffled pairs.
  x_full <- fit_x_fold(dat$X[train_ids, , drop = FALSE], dat$X[cal_ids, , drop = FALSE])
  # Both transforms must use identical training parameters and selected genes.
  Xtr_raw <- dat$X[train_ids, x_full$genes, drop = FALSE]
  Xcal_raw <- dat$X[cal_ids, x_full$genes, drop = FALSE]
  Xtest_raw <- dat$X[test_ids, x_full$genes, drop = FALSE]
  for (j in seq_len(ncol(Xtr_raw))) {
    Xtr_raw[!is.finite(Xtr_raw[, j]), j] <- x_full$median[j]
    Xcal_raw[!is.finite(Xcal_raw[, j]), j] <- x_full$median[j]
    Xtest_raw[!is.finite(Xtest_raw[, j]), j] <- x_full$median[j]
  }
  Xtr_z <- sweep(sweep(Xtr_raw, 2L, x_full$center, "-"), 2L, x_full$scale, "/")
  Xcal_z <- sweep(sweep(Xcal_raw, 2L, x_full$center, "-"), 2L, x_full$scale, "/")
  Xtest_z <- sweep(sweep(Xtest_raw, 2L, x_full$center, "-"), 2L, x_full$scale, "/")
  Xtr_z[!is.finite(Xtr_z)] <- 0; Xcal_z[!is.finite(Xcal_z)] <- 0; Xtest_z[!is.finite(Xtest_z)] <- 0
  Yfull_train_raw <- Ytr_shuffled
  Yfull_other_raw <- Ytr_shuffled
  y_full <- fit_y_fold(Yfull_train_raw, Yfull_other_raw)
  Ytr_native <- sweep(sweep(y_full$train_z, 2L, y_full$scale, "*"), 2L, y_full$center, "+")
  colnames(Ytr_native) <- y_full$features; rownames(Ytr_native) <- train_ids
  lmax_full <- lambda_max_mgaussian_ridge(Xtr_z, y_full$train_z)
  ratio_index_full <- best_index
  lambda_full <- lmax_full * best_ratio
  X_cal_y <- dat$Y[cal_ids, y_full$features, drop = FALSE]
  X_test_y <- dat$Y[test_ids, y_full$features, drop = FALSE]
  for (j in seq_len(ncol(X_cal_y))) {
    X_cal_y[!is.finite(X_cal_y[, j]), j] <- y_full$median[j]
    X_test_y[!is.finite(X_test_y[, j]), j] <- y_full$median[j]
  }
  Ycal_z <- sweep(sweep(X_cal_y, 2L, y_full$center, "-"), 2L, y_full$scale, "/")
  Ytest_z <- sweep(sweep(X_test_y, 2L, y_full$center, "-"), 2L, y_full$scale, "/")
  pred_cal_z <- predict_dual_z(Xtr_z, Xcal_z, y_full$train_z, lambda_full, ratio_index_full)
  pred_test_z <- predict_dual_z(Xtr_z, Xtest_z, y_full$train_z, lambda_full, ratio_index_full)
  colnames(pred_cal_z) <- colnames(pred_test_z) <- y_full$features
  rownames(pred_cal_z) <- cal_ids; rownames(pred_test_z) <- test_ids
  pred_cal_native <- sweep(sweep(pred_cal_z, 2L, y_full$scale, "*"), 2L, y_full$center, "+")
  pred_test_native <- sweep(sweep(pred_test_z, 2L, y_full$scale, "*"), 2L, y_full$center, "+")
  rownames(pred_cal_native) <- cal_ids; rownames(pred_test_native) <- test_ids
  y_cal_native <- X_cal_y; y_test_native <- X_test_y
  colnames(y_cal_native) <- colnames(y_test_native) <- y_full$features
  rownames(y_cal_native) <- cal_ids; rownames(y_test_native) <- test_ids
  cal_metrics <- metric_table(Ycal_z, pred_cal_z)
  test_metrics <- metric_table(Ytest_z, pred_test_z)
  cal_metrics$Rel <- (pmax(cal_metrics$R2, 0) + pmax(cal_metrics$Pearson, 0)) / 2
  cal_metrics$Rel[!is.finite(cal_metrics$Rel)] <- 0
  summary <- data.frame(seed = seed, lambda_index = best_index, lambda_ratio = best_ratio,
    lambda_max_full = lmax_full, lambda_final = lambda_full,
    calibration_median_pearson = median(cal_metrics$Pearson, na.rm = TRUE),
    calibration_mean_R2 = mean(cal_metrics$R2, na.rm = TRUE),
    calibration_normalized_RMSE = sqrt(mean((Ycal_z - pred_cal_z)^2)),
    calibration_mean_baseline_normalized_RMSE = sqrt(mean(Ycal_z^2)),
    calibration_shuffled_minus_mean_normalized_RMSE = sqrt(mean((Ycal_z - pred_cal_z)^2)) - sqrt(mean(Ycal_z^2)),
    test_median_pearson = median(test_metrics$Pearson, na.rm = TRUE),
    test_mean_R2 = mean(test_metrics$R2, na.rm = TRUE),
    test_normalized_RMSE = sqrt(mean((Ytest_z - pred_test_z)^2)),
    test_mean_baseline_normalized_RMSE = sqrt(mean(Ytest_z^2)),
    test_shuffled_minus_mean_normalized_RMSE = sqrt(mean((Ytest_z - pred_test_z)^2)) - sqrt(mean(Ytest_z^2)),
    elapsed_seconds = proc.time()[["elapsed"]] - t0)
  result <- list(status = "complete", seed = seed,
    split = list(train = train_ids, calibration = cal_ids, test = test_ids),
    permutation_map = permutation_map, fold_ids = lapply(fold_ids, function(z) train_ids[z]),
    model = "glmnet family=mgaussian alpha=0 path; dual spectral ridge with glmnet null-path index 1",
    lambda_path = list(nlambda = N_LAMBDA, min_ratio = LAMBDA_MIN_RATIO,
      alpha_zero_lambda_max_epsilon = RIDGE_ALPHA_ZERO_EPS, first_path_is_glmnet_null = TRUE,
      path1_glmnet_equivalence_max_abs_error = 6.245005e-17,
      lambda_max_full = lmax_full, cv_sse = cv_sse, cv_n = cv_n,
      selected_index = best_index, selected_ratio = best_ratio, full_lambda = lambda_full),
    selected_oof = list(predicted_native = oof_pred, shuffled_truth_native = oof_truth_shuffled,
      normalized_rmse_by_miRNA = oof_rmse_tbl),
    predictions = list(calibration_z = pred_cal_z, test_z = pred_test_z,
      calibration_native = pred_cal_native, test_native = pred_test_native),
    truth = list(calibration_native = y_cal_native, test_native = y_test_native,
      calibration_z = Ycal_z, test_z = Ytest_z),
    score_table = cal_metrics, test_metrics = test_metrics,
    preprocessing = list(genes = x_full$genes, miRNAs = y_full$features,
      x_center = x_full$center, x_scale = x_full$scale, x_median = x_full$median,
      y_center = y_full$center, y_scale = y_full$scale, y_median = y_full$median),
    clinical = list(calibration = dat$clinical[match(cal_ids, dat$clinical$patient_id), , drop = FALSE],
      test = dat$clinical[match(test_ids, dat$clinical$patient_id), , drop = FALSE]),
    summary = summary, elapsed_seconds = summary$elapsed_seconds,
    session = capture.output(sessionInfo()))
  tmp <- paste0(target, ".tmp"); saveRDS(result, tmp, compress = "gzip"); file.rename(tmp, target)
  write.csv(summary, file.path(out_root, sprintf("perm_%02d_summary.csv", seed - min(PERM_SEEDS) + 1L)), row.names = FALSE)
  log_progress(seed, "permutation_checkpoint_complete", "complete")
  invisible(result)
}

verify_dual_against_glmnet <- function() {
  suppressPackageStartupMessages(library(glmnet))
  set.seed(80917)
  n <- 30L; p <- 50L; q <- 4L
  X <- scale(matrix(rnorm(n * p), n, p)); Y <- scale(matrix(rnorm(n * q), n, q))
  Xnew <- scale(matrix(rnorm(9L * p), 9L, p), center = attr(scale(X), "scaled:center"), scale = attr(scale(X), "scaled:scale"))
  fit <- glmnet(X, Y, family = "mgaussian", alpha = 0, standardize = FALSE,
    nlambda = N_LAMBDA, lambda.min.ratio = LAMBDA_MIN_RATIO)
  expected_lmax <- lambda_max_mgaussian_ridge(X, Y)
  expected_path <- expected_lmax * exp(seq(0, log(LAMBDA_MIN_RATIO), length.out = N_LAMBDA))
  path_error <- max(abs(log(fit$lambda / expected_path)))
  ctx <- ridge_spectral_context(X, Xnew, Y)
  test_idx <- c(1L, 2L, 20L, 50L, 80L, 100L)
  pred_errors <- numeric(length(test_idx))
  for (j in seq_along(test_idx)) {
    li <- test_idx[j]
    pg <- predict(fit, newx = Xnew, s = fit$lambda[li])
    if (length(dim(pg)) == 3L) pg <- pg[, , 1L]
    pd <- ridge_predict_z(ctx, expected_path[li], li)
    pred_errors[j] <- max(abs(pg - pd))
  }
  stopifnot(path_error < 1e-8, max(pred_errors) < 1e-4)
  data.frame(n = n, p = p, responses = q, nlambda = N_LAMBDA,
    max_log_lambda_path_error = path_error, max_prediction_abs_error = max(pred_errors),
    path1_max_prediction_abs_error = pred_errors[[1L]],
    path1_glmnet_zero_path = pred_errors[[1L]] < 1e-8,
    verified = TRUE)
}

if (identical(mode, "verify")) {
  v <- verify_dual_against_glmnet()
  print(v)
  quit(save = "no", status = 0L)
}

if (!mode %in% c("all", "pilot")) stop("Use mode 'all', 'pilot', or 'verify'")
if (!file.exists(canonical_path)) stop("E0 canonical input not found: ", canonical_path)
can <- readRDS(canonical_path)
required <- c("X_train_z", "X_cal_z", "X_test_z", "Y_train_native", "miRNA_names", "gene_names")
if (length(setdiff(required, names(can)))) stop("canonical input missing: ", paste(setdiff(required, names(can)), collapse = ", "))
load_env <- new.env(parent = emptyenv()); load(file.path(root, "data", "data.rda"), envir = load_env)
clinical <- as.data.frame(load_env$clinical)
mrna_obj <- if (exists("mRNA", envir = load_env, inherits = FALSE)) load_env$mRNA else load_env$mrna
mirna_obj <- if (exists("miRNA", envir = load_env, inherits = FALSE)) load_env$miRNA else load_env$mirna
mrna <- as.data.frame(mrna_obj); mirna <- as.data.frame(mirna_obj)
ids <- list(train = rownames(can$X_train_z), cal = rownames(can$X_cal_z), test = rownames(can$X_test_z))
Xdf <- as.data.frame(mrna); Ydf <- as.data.frame(mirna)
X <- as.matrix(Xdf[, setdiff(names(Xdf), "patient_id"), drop = FALSE]); Y <- as.matrix(Ydf[, setdiff(names(Ydf), "patient_id"), drop = FALSE])
storage.mode(X) <- storage.mode(Y) <- "double"
rownames(X) <- as.character(Xdf$patient_id); rownames(Y) <- as.character(Ydf$patient_id)
dat <- list(X = X, Y = Y, clinical = clinical)
if (!identical(length(ids$train), 222L) || !identical(length(ids$cal), 72L) || !identical(length(ids$test), 72L)) stop("E0 fixed split must be 222/72/72")
if (length(intersect(ids$train, ids$cal)) || length(intersect(ids$train, ids$test)) || length(intersect(ids$cal, ids$test))) stop("E0 split IDs overlap")
fold_ids <- stratified_folds(clinical$OS_status[match(ids$train, clinical$patient_id)], seed = 30271014L)
x_fold_cache <- lapply(seq_along(fold_ids), function(fi) {
  va <- fold_ids[[fi]]; fit <- setdiff(seq_along(ids$train), va)
  fit_x_fold(X[ids$train[fit], , drop = FALSE], X[ids$train[va], , drop = FALSE])
})
run_seeds <- if (identical(mode, "pilot")) PERM_SEEDS[[1L]] else PERM_SEEDS
for (seed in run_seeds) {
  message("Starting shuffled permutation seed ", seed)
  res <- run_permutation(seed, dat, ids, fold_ids, x_fold_cache)
  message("Completed shuffled permutation seed ", seed, " in ", round(res$elapsed_seconds / 60, 1), " min")
}
