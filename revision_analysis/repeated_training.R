#!/usr/bin/env Rscript
# E1 repeated patient-level training and evaluation.
# This script starts from data/data.rda on every run and never reads the fixed-split
# matrices in output/gpu_model_inputs or new/results/E0.

suppressPackageStartupMessages(library(glmnet))

args <- commandArgs(trailingOnly = TRUE)
seed_arg <- if (length(args) >= 1L) args[[1L]] else "all"
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_path <- file.path(root, "data", "data.rda")
out_root <- file.path(root, "new", "results", "E1")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
progress_path <- file.path(out_root, "current_progress.csv")

log_progress <- function(seed, model = NA_character_, config = NA_integer_, foldstage, status = "running", elapsed_seconds = NA_real_) {
  row <- data.frame(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), seed = as.integer(seed),
    model = as.character(model), config = as.integer(config), foldstage = as.character(foldstage),
    status = as.character(status), elapsed_seconds = as.numeric(elapsed_seconds), check.names = FALSE)
  write.table(row, progress_path, sep = ",", row.names = FALSE, col.names = !file.exists(progress_path),
    append = file.exists(progress_path), na = "NA", qmethod = "double")
  invisible(NULL)
}

seeds <- 20261014:20261033
if (!identical(seed_arg, "all") && !identical(seed_arg, "sanity") && !identical(seed_arg, "integration-sanity")) {
  requested <- as.integer(strsplit(seed_arg, ",", fixed = TRUE)[[1L]])
  if (anyNA(requested) || any(!requested %in% seeds)) stop("seed argument must be 'all' or comma-separated seeds 20261014..20261033")
  seeds <- requested
}

K_GRID <- c(20L, 50L, 100L, 150L, 200L)
ALPHA_GRID <- c(0.1, 0.3, 0.5, 0.7, 0.9)
TOP_MRNA <- 3000L
N_FOLDS <- 5L
PROTOCOL_VERSION <- "trainCV_lambda_default_v2"
N_LAMBDA <- 100L
LAMBDA_MIN_RATIO_POLICY <- "0.01 when n < p; 0.0001 when n >= p (glmnet default)"

load_inputs <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  required <- c("clinical", "mrna", "mirna")
  if (!all(required %in% ls(env))) stop("data.rda must contain clinical, mrna, and mirna")
  clinical <- as.data.frame(env$clinical)
  x_df <- as.data.frame(env$mrna)
  y_df <- as.data.frame(env$mirna)
  id <- "patient_id"
  if (!all(vapply(list(clinical, x_df, y_df), function(d) id %in% names(d), logical(1)))) stop("patient_id missing from an input table")
  common <- Reduce(intersect, list(as.character(clinical[[id]]), as.character(x_df[[id]]), as.character(y_df[[id]])))
  if (length(common) != 366L) warning("Matched patient count is ", length(common), "; expected 366")
  clinical <- clinical[match(common, as.character(clinical[[id]])), , drop = FALSE]
  x_df <- x_df[match(common, as.character(x_df[[id]])), , drop = FALSE]
  y_df <- y_df[match(common, as.character(y_df[[id]])), , drop = FALSE]
  X <- as.matrix(x_df[, setdiff(names(x_df), id), drop = FALSE])
  Y <- as.matrix(y_df[, setdiff(names(y_df), id), drop = FALSE])
  storage.mode(X) <- "double"; storage.mode(Y) <- "double"
  rownames(X) <- rownames(Y) <- common
  list(clinical = clinical, X = X, Y = Y, ids = common)
}

stratified_partition <- function(group, sizes, seed) {
  stopifnot(sum(sizes) == length(group))
  set.seed(seed)
  group <- as.character(group); group[is.na(group) | !nzchar(group)] <- "<NA>"
  pieces <- vector("list", length(sizes)); names(pieces) <- names(sizes)
  for (g in unique(group)) {
    ix <- sample(which(group == g)); n <- length(ix)
    exact <- n * sizes / sum(sizes)
    alloc <- floor(exact)
    rem <- n - sum(alloc)
    if (rem > 0L) alloc[order(exact - alloc, decreasing = TRUE)[seq_len(rem)]] <- alloc[order(exact - alloc, decreasing = TRUE)[seq_len(rem)]] + 1L
    start <- 1L
    for (j in seq_along(sizes)) {
      if (alloc[j] > 0L) pieces[[j]] <- c(pieces[[j]], ix[start:(start + alloc[j] - 1L)])
      start <- start + alloc[j]
    }
  }
  lapply(pieces, sort)
}

make_folds <- function(group, k, seed) {
  set.seed(seed)
  group <- as.character(group); group[is.na(group) | !nzchar(group)] <- "<NA>"
  folds <- vector("list", k)
  for (g in unique(group)) {
    ix <- sample(which(group == g)); assignment <- rep(seq_len(k), length.out = length(ix))
    for (f in seq_len(k)) folds[[f]] <- c(folds[[f]], ix[assignment == f])
  }
  lapply(folds, sort)
}

fit_preprocess <- function(X, Y, top_mrna = TOP_MRNA) {
  finite_x <- is.finite(X); finite_y <- is.finite(Y)
  keep_x <- colMeans(finite_x) >= 0.8; keep_y <- colMeans(finite_y) >= 0.8
  X <- X[, keep_x, drop = FALSE]; Y <- Y[, keep_y, drop = FALSE]
  med_x <- apply(X, 2L, function(v) { v <- v[is.finite(v)]; if (!length(v)) 0 else median(v) })
  med_y <- apply(Y, 2L, function(v) { v <- v[is.finite(v)]; if (!length(v)) 0 else median(v) })
  for (j in seq_len(ncol(X))) X[!is.finite(X[, j]), j] <- med_x[j]
  for (j in seq_len(ncol(Y))) Y[!is.finite(Y[, j]), j] <- med_y[j]
  xv <- apply(X, 2L, var); yv <- apply(Y, 2L, var)
  kx <- is.finite(xv) & xv > 1e-8; ky <- is.finite(yv) & yv > 1e-8
  X <- X[, kx, drop = FALSE]; Y <- Y[, ky, drop = FALSE]
  xv <- xv[kx]
  genes <- names(sort(xv, decreasing = TRUE))[seq_len(min(top_mrna, length(xv)))]
  X <- X[, genes, drop = FALSE]
  x_center <- colMeans(X); x_scale <- apply(X, 2L, sd); x_scale[!is.finite(x_scale) | x_scale < 1e-8] <- 1
  y_center <- colMeans(Y); y_scale <- apply(Y, 2L, sd); y_scale[!is.finite(y_scale) | y_scale < 1e-8] <- 1
  list(genes = genes, mirnas = colnames(Y), x_median = med_x[genes], y_median = med_y[colnames(Y)],
       x_center = x_center, x_scale = x_scale, y_center = y_center, y_scale = y_scale)
}

apply_preprocess <- function(pp, X, Y = NULL) {
  X <- X[, pp$genes, drop = FALSE]
  for (j in seq_len(ncol(X))) X[!is.finite(X[, j]), j] <- pp$x_median[j]
  Xz <- sweep(sweep(X, 2L, pp$x_center, "-"), 2L, pp$x_scale, "/"); Xz[!is.finite(Xz)] <- 0
  if (is.null(Y)) return(list(X = Xz))
  Y <- Y[, pp$mirnas, drop = FALSE]
  for (j in seq_len(ncol(Y))) Y[!is.finite(Y[, j]), j] <- pp$y_median[j]
  Yz <- sweep(sweep(Y, 2L, pp$y_center, "-"), 2L, pp$y_scale, "/"); Yz[!is.finite(Yz)] <- 0
  list(X = Xz, Y = Yz, Y_native = Y)
}

metric_summary <- function(real, pred, scale = rep(1, ncol(real))) {
  r <- vapply(seq_len(ncol(real)), function(j) suppressWarnings(cor(real[, j], pred[, j], use = "complete.obs")), numeric(1))
  denom <- colSums(sweep(real, 2L, colMeans(real), "-")^2)
  sse <- colSums((real - pred)^2)
  r2 <- 1 - sse / denom
  rmse <- sqrt(colMeans((real - pred)^2))
  normalized_error <- if (is.matrix(scale)) {
    (real - pred) / scale
  } else {
    sweep(real - pred, 2L, scale, "/")
  }
  c(median_pearson = median(r[is.finite(r)], na.rm = TRUE), mean_r2 = mean(r2[is.finite(r2)], na.rm = TRUE),
    overall_rmse = sqrt(mean((real - pred)^2)), normalized_rmse = sqrt(mean(normalized_error^2)),
    median_mirna_rmse = median(rmse[is.finite(rmse)], na.rm = TRUE),
    n_valid_cor = sum(is.finite(r)), n_valid_r2 = sum(is.finite(r2)), negative_r2_fraction = mean(r2 < 0, na.rm = TRUE))
}

fit_glmnet_path <- function(x, y, alpha) {
  lambda_min_ratio <- if (nrow(x) < ncol(x)) 0.01 else 0.0001
  glmnet(x = x, y = y, family = "mgaussian", alpha = alpha, standardize = FALSE,
         nlambda = N_LAMBDA, lambda.min.ratio = lambda_min_ratio)
}

predict_native <- function(fit, x, lambda, pp, row_ids) {
  p <- predict(fit, newx = x, s = lambda)
  if (length(dim(p)) == 3L) p <- p[, , 1L]
  p <- as.matrix(p); colnames(p) <- pp$mirnas; rownames(p) <- row_ids
  sweep(sweep(p, 2L, pp$y_scale, "*"), 2L, pp$y_center, "+")
}

run_seed <- function(seed, dat) {
  target <- file.path(out_root, sprintf("seed_%d.rds", seed))
  if (file.exists(target)) {
    old <- tryCatch(readRDS(target), error = function(e) NULL)
    canonical_ok <- !is.null(old$canonical_paths) && all(c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge") %in% names(old$canonical_paths)) &&
      all(file.exists(old$canonical_paths[c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")]))
    if (!is.null(old) && identical(old$seed, seed) && identical(old$status, "complete") &&
        identical(old$protocol_version, PROTOCOL_VERSION) && canonical_ok) {
      log_progress(seed, foldstage = "seed_complete_artifacts_verified", status = "skipped_complete", elapsed_seconds = old$elapsed_seconds)
      return(invisible(old))
    }
  }
  t0 <- proc.time()[["elapsed"]]
  log_progress(seed, foldstage = "seed_started")
  set.seed(seed)
  partitions <- stratified_partition(dat$clinical$OS_status, c(train = 222L, calibration = 72L, test = 72L), seed)
  tr <- partitions$train; ca <- partitions$calibration; te <- partitions$test
  folds <- make_folds(dat$clinical$OS_status[tr], N_FOLDS, seed + 10000L)
  configs <- list(ElasticNet = expand.grid(alpha = ALPHA_GRID), PCA_ElasticNet = expand.grid(k = K_GRID, alpha = ALPHA_GRID), PCA_Ridge = data.frame(k = K_GRID, alpha = 0))
  # Per-model/config OOF predictions are restored to the original response scale
  # before pooling across folds. Lambda is selected only by training OOF RMSE.
  # Preprocessing and PCA are fitted once per fold and reused across all model
  # configurations; this avoids repeating the expensive 18k-feature transforms.
  fold_cache <- vector("list", length(folds))
  for (fi in seq_along(folds)) {
    va_local <- folds[[fi]]; fit_local <- setdiff(seq_along(tr), va_local)
    id_fit <- tr[fit_local]; id_val <- tr[va_local]
    pp <- fit_preprocess(dat$X[id_fit, , drop = FALSE], dat$Y[id_fit, , drop = FALSE])
    a <- apply_preprocess(pp, dat$X[id_fit, , drop = FALSE], dat$Y[id_fit, , drop = FALSE])
    b <- apply_preprocess(pp, dat$X[id_val, , drop = FALSE], dat$Y[id_val, , drop = FALSE])
    fold_cache[[fi]] <- list(rows = va_local, pp = pp, train = a, valid = b,
      pca = prcomp(a$X, center = FALSE, scale. = FALSE))
    log_progress(seed, foldstage = paste0("preprocess_pca_fold_", fi, "_of_", N_FOLDS), elapsed_seconds = proc.time()[["elapsed"]] - t0)
  }
  min_fold_pca_rank <- min(vapply(fold_cache, function(fc) ncol(fc$pca$x), integer(1)))
  k_rank_applicability <- data.frame(
    model = rep(c("PCA_ElasticNet", "PCA_Ridge"), each = length(K_GRID)),
    k = rep(K_GRID, times = 2L),
    minimum_available_fold_rank = min_fold_pca_rank,
    applicable = rep(K_GRID <= min_fold_pca_rank, times = 2L)
  )
  score_rows <- list(); oof_cache <- list()
  for (model in names(configs)) {
    cfg <- configs[[model]]
    for (ci in seq_len(nrow(cfg))) {
      log_progress(seed, model, ci, "config_started", elapsed_seconds = proc.time()[["elapsed"]] - t0)
      vals <- vector("list", length(folds)); applicable <- TRUE
      for (fi in seq_along(folds)) {
        fc <- fold_cache[[fi]]; va_local <- fc$rows; pp <- fc$pp; a <- fc$train; b <- fc$valid
        if (model == "ElasticNet") {
          x_fit <- a$X; x_val <- b$X; k <- NA_integer_; alpha <- cfg$alpha[ci]
        } else {
          pca <- fc$pca
          k <- cfg$k[ci]
          if (k >= min(dim(a$X)) || k > ncol(pca$x)) { applicable <- FALSE; break }
          x_fit <- pca$x[, seq_len(k), drop = FALSE]
          x_val <- predict(pca, newdata = b$X)[, seq_len(k), drop = FALSE]
          alpha <- cfg$alpha[ci]
        }
        fit <- tryCatch(fit_glmnet_path(x_fit, a$Y, alpha), error = function(e) NULL)
        if (is.null(fit)) { applicable <- FALSE; break }
        lambdas <- fit$lambda
        p <- predict(fit, newx = x_val, s = lambdas)
        if (length(dim(p)) == 3L) p <- aperm(p, c(1L, 2L, 3L))
        # predict(mgaussian) is samples x responses x lambdas.
        native <- array(NA_real_, dim = c(nrow(x_val), ncol(a$Y), length(lambdas)))
        for (li in seq_along(lambdas)) {
          pl <- if (length(dim(p)) == 3L) p[, , li] else as.matrix(p)
          native[, , li] <- sweep(sweep(pl, 2L, pp$y_scale, "*"), 2L, pp$y_center, "+")
        }
        vals[[fi]] <- list(rows = va_local, real = b$Y_native, pred = native, lambdas = lambdas,
                           lambda_min_ratio = if (nrow(x_fit) < ncol(x_fit)) 0.01 else 0.0001,
                           lambda_min_ratio_observed = tail(lambdas, 1L) / lambdas[1L],
                           nlambda_observed = length(lambdas),
                           y_scale = pp$y_scale, names = pp$mirnas)
        log_progress(seed, model, ci, paste0("fold_", fi, "_of_", N_FOLDS, "_glmnet_complete"), elapsed_seconds = proc.time()[["elapsed"]] - t0)
        rm(fit, p, native)
      }
      if (!applicable) {
        log_progress(seed, model, ci, "config_inapplicable_or_fit_failed", status = "skipped", elapsed_seconds = proc.time()[["elapsed"]] - t0)
        rm(vals); invisible(gc(FALSE))
        next
      }
      common_y <- Reduce(intersect, lapply(vals, function(v) v$names))
      if (length(common_y) < 2L) next
      # Select a common relative lambda-path position. Absolute lambda scales
      # differ by fold and cannot be averaged then transferred to a full refit.
      n_path <- min(vapply(vals, function(v) length(v$lambdas), integer(1)))
      errors <- rep(0, n_path); counts <- rep(0, n_path)
      for (v in vals) {
        yi <- match(common_y, v$names)
        for (li in seq_len(n_path)) {
          residual <- v$real[, yi, drop = FALSE] - v$pred[, yi, li, drop = FALSE][, , 1L]
          residual_z <- sweep(residual, 2L, v$y_scale[yi], "/")
          errors[li] <- errors[li] + sum(residual_z^2)
          counts[li] <- counts[li] + length(v$rows) * length(yi)
        }
      }
      best_path_index <- which.min(errors / counts)
      picked_ratio <- vapply(vals, function(v) v$lambdas[best_path_index] / v$lambdas[1L], numeric(1))
      chosen_ratio <- exp(mean(log(picked_ratio)))
      oof_real <- matrix(NA_real_, length(tr), length(common_y), dimnames = list(dat$ids[tr], common_y))
      oof_pred <- oof_real
      oof_scale <- oof_real
      for (fi in seq_along(vals)) {
        v <- vals[[fi]]; li <- best_path_index
        yi <- match(common_y, v$names)
        oof_real[v$rows, ] <- v$real[, yi, drop = FALSE]
        oof_pred[v$rows, ] <- v$pred[, yi, li, drop = FALSE][, , 1L]
        oof_scale[v$rows, ] <- matrix(rep(v$y_scale[yi], each = length(v$rows)), nrow = length(v$rows))
      }
      key <- paste(model, ci, sep = "__")
      oof_norm_rmse <- sqrt(colMeans(((oof_real - oof_pred) / oof_scale)^2))
      oof_cache[[key]] <- list(real = oof_real, pred = oof_pred, scale = oof_scale,
                               normalized_rmse_by_miRNA = oof_norm_rmse,
                               lambda_path_index = best_path_index, lambda_ratio = chosen_ratio,
                               config = cfg[ci, , drop = FALSE])
      ms <- metric_summary(oof_real, oof_pred, scale = oof_scale)
      score_rows[[length(score_rows) + 1L]] <- data.frame(model = model, config = ci, k = if ("k" %in% names(cfg)) cfg$k[ci] else NA_integer_,
        alpha = cfg$alpha[ci], lambda_path_index = best_path_index, lambda_ratio = chosen_ratio,
        nlambda_requested = N_LAMBDA, cv_common_nlambda = n_path,
        lambda_min_ratio_requested_by_shape = paste(sort(unique(vapply(vals, `[[`, numeric(1), "lambda_min_ratio"))), collapse = ";"),
        lambda_min_ratio_observed_by_fold = paste(format(vapply(vals, `[[`, numeric(1), "lambda_min_ratio_observed"), scientific = TRUE), collapse = ";"),
        nlambda_observed_by_fold = paste(vapply(vals, `[[`, integer(1), "nlambda_observed"), collapse = ";"),
        lambda_selection_rmse_scale = "fold_training_SD_normalized", t(ms), check.names = FALSE)
      log_progress(seed, model, ci, "config_oof_scored", status = "complete", elapsed_seconds = proc.time()[["elapsed"]] - t0)
      rm(vals, oof_real, oof_pred, oof_scale, errors, counts)
      invisible(gc(FALSE))
    }
  }
  tuning <- if (length(score_rows)) do.call(rbind, score_rows) else stop("No valid OOF configurations")
  # Apply the pre-registered lexicographic rule: median Pearson, mean R2, normalized RMSE.
  selected <- lapply(base::split(tuning, tuning$model), function(d) {
    d <- d[order(-d$median_pearson, -d$mean_r2, d$normalized_rmse), , drop = FALSE]; d[1L, , drop = FALSE]
  })
  model_order <- names(selected)
  aggregate <- do.call(rbind, lapply(selected, function(d) d[, c("model", "config", "k", "alpha", "lambda_path_index", "lambda_ratio", "median_pearson", "mean_r2", "normalized_rmse", "overall_rmse"), drop = FALSE]))
  aggregate <- aggregate[order(-aggregate$median_pearson, -aggregate$mean_r2, aggregate$normalized_rmse), , drop = FALSE]
  winner <- aggregate$model[1L]
  selected_oof <- lapply(selected, function(row) {
    key <- paste(row$model, row$config, sep = "__")
    oof_cache[[key]]
  })
  names(selected_oof) <- names(selected)
  oof_cache <- NULL
  invisible(gc(FALSE))
  # Refit each selected model on all 222 training cases; calibration and test are prediction-only.
  pp <- fit_preprocess(dat$X[tr, , drop = FALSE], dat$Y[tr, , drop = FALSE])
  trp <- apply_preprocess(pp, dat$X[tr, , drop = FALSE], dat$Y[tr, , drop = FALSE])
  cap <- apply_preprocess(pp, dat$X[ca, , drop = FALSE], dat$Y[ca, , drop = FALSE])
  tep <- apply_preprocess(pp, dat$X[te, , drop = FALSE], dat$Y[te, , drop = FALSE])
  pred <- list(calibration = list(), test = list()); model_meta <- list()
  for (model in names(selected)) {
    row <- selected[[model]]; key <- paste(model, row$config, sep = "__"); lambda_ratio <- row$lambda_ratio
    if (model == "ElasticNet") {
      xx <- trp$X; xc <- cap$X; xt <- tep$X
    } else {
      pca <- prcomp(trp$X, center = FALSE, scale. = FALSE); k <- as.integer(row$k)
      xx <- pca$x[, seq_len(k), drop = FALSE]
      xc <- predict(pca, newdata = cap$X)[, seq_len(k), drop = FALSE]
      xt <- predict(pca, newdata = tep$X)[, seq_len(k), drop = FALSE]
    }
    alpha <- row$alpha
    fit <- fit_glmnet_path(xx, trp$Y, alpha)
    final_path_index <- which.min(abs(log(fit$lambda / fit$lambda[1L]) - log(lambda_ratio)))
    lambda <- fit$lambda[final_path_index]
    coef_parts <- coef(fit, s = lambda)
    coef_matrix <- do.call(cbind, lapply(coef_parts, as.matrix))
    colnames(coef_matrix) <- pp$mirnas
    pred$calibration[[model]] <- predict_native(fit, xc, lambda, pp, dat$ids[ca])
    pred$test[[model]] <- predict_native(fit, xt, lambda, pp, dat$ids[te])
    model_meta[[model]] <- list(config = row, lambda = lambda, lambda_path_index = final_path_index,
      lambda_ratio_requested = lambda_ratio, lambda_ratio_final = fit$lambda[final_path_index] / fit$lambda[1L],
      lambda_path = list(nlambda_requested = N_LAMBDA, nlambda_observed = length(fit$lambda),
        lambda_min_ratio_requested = if (nrow(xx) < ncol(xx)) 0.01 else 0.0001,
        lambda_min_ratio_observed = tail(fit$lambda, 1L) / fit$lambda[1L],
        policy = LAMBDA_MIN_RATIO_POLICY, selection = "training-only five-fold OOF relative lambda path"),
      pca_k = if (model == "ElasticNet") NA_integer_ else as.integer(row$k), coefficients = coef_matrix,
      pca = if (model == "ElasticNet") NULL else pca)
    log_progress(seed, model, row$config, "final_refit_and_cal_test_prediction_complete", status = "complete", elapsed_seconds = proc.time()[["elapsed"]] - t0)
  }
  # Mean baseline is trained only from the 222 training responses.
  pred$calibration$Mean <- matrix(rep(pp$y_center, each = length(ca)), nrow = length(ca), dimnames = list(dat$ids[ca], pp$mirnas))
  pred$test$Mean <- matrix(rep(pp$y_center, each = length(te)), nrow = length(te), dimnames = list(dat$ids[te], pp$mirnas))
  # Calibration-only reliability for the training-selected winner.
  cal_real <- cap$Y_native; cal_hat <- pred$calibration[[winner]][, colnames(cal_real), drop = FALSE]
  cal_r <- vapply(seq_len(ncol(cal_real)), function(j) suppressWarnings(cor(cal_real[, j], cal_hat[, j], use = "complete.obs")), numeric(1))
  cal_sse <- colSums((cal_real - cal_hat)^2); cal_den <- colSums(sweep(cal_real, 2L, colMeans(cal_real), "-")^2)
  cal_r2 <- 1 - cal_sse / cal_den
  reliability <- data.frame(miRNA = colnames(cal_real), Pearson = cal_r, R2 = cal_r2,
    Rel = (pmax(cal_r, 0) + pmax(cal_r2, 0)) / 2, RMSE = sqrt(colMeans((cal_real - cal_hat)^2)))
  test_metrics <- do.call(rbind, lapply(names(pred$test), function(m) {
    z <- metric_summary(tep$Y_native[, colnames(pred$test[[m]]), drop = FALSE], pred$test[[m]], scale = pp$y_scale[colnames(pred$test[[m]])])
    data.frame(model = m, t(z), check.names = FALSE)
  }))
  # Per-model OOF diagnostics are training-only and remain separated from
  # calibration/test evaluation so they can serve as E2/E3 OOF baselines.
  selected_oof <- lapply(selected_oof, function(z) {
    z$normalized_rmse_by_miRNA <- data.frame(miRNA = names(z$normalized_rmse_by_miRNA),
      OOF_RMSE_normalized = as.numeric(z$normalized_rmse_by_miRNA))
    z
  })
  canonical_paths <- character()
  for (model in names(pred$calibration)[names(pred$calibration) != "Mean"]) {
    y_cal_z <- sweep(sweep(cap$Y_native, 2L, pp$y_center, "-"), 2L, pp$y_scale, "/")
    y_test_z <- sweep(sweep(tep$Y_native, 2L, pp$y_center, "-"), 2L, pp$y_scale, "/")
    pred_cal_z <- sweep(sweep(pred$calibration[[model]], 2L, pp$y_center, "-"), 2L, pp$y_scale, "/")
    pred_test_z <- sweep(sweep(pred$test[[model]], 2L, pp$y_center, "-"), 2L, pp$y_scale, "/")
    r <- vapply(seq_len(ncol(y_cal_z)), function(j) suppressWarnings(cor(y_cal_z[, j], pred_cal_z[, j], use = "complete.obs")), numeric(1))
    r2 <- 1 - colSums((y_cal_z - pred_cal_z)^2) / colSums(sweep(y_cal_z, 2L, colMeans(y_cal_z), "-")^2)
    score_table <- data.frame(miRNA = colnames(y_cal_z), R2 = r2, Pearson = r,
      RMSE = sqrt(colMeans((y_cal_z - pred_cal_z)^2)))
    score_table$Rel <- (pmax(score_table$R2, 0) + pmax(score_table$Pearson, 0)) / 2
    score_table$Rel[!is.finite(score_table$Rel)] <- 0
    canonical <- list(protocol_version = PROTOCOL_VERSION,
      lambda_path_policy = list(nlambda = N_LAMBDA, min_ratio_by_shape = LAMBDA_MIN_RATIO_POLICY),
      X_train_z = trp$X, X_cal_z = cap$X, X_test_z = tep$X,
      Y_train_z = trp$Y, Y_cal_z = y_cal_z, Y_test_z = y_test_z,
      X_train_native = dat$X[tr, pp$genes, drop = FALSE], X_cal_native = dat$X[ca, pp$genes, drop = FALSE],
      X_test_native = dat$X[te, pp$genes, drop = FALSE],
      Y_train_native = trp$Y_native, Y_cal_native = cap$Y_native, Y_test_native = tep$Y_native,
      pred_cal_z = pred_cal_z, pred_test_z = pred_test_z,
      pred_cal_native = pred$calibration[[model]], pred_test_native = pred$test[[model]],
      clinical_train = dat$clinical[tr, , drop = FALSE], clinical_cal = dat$clinical[ca, , drop = FALSE],
      clinical_test = dat$clinical[te, , drop = FALSE],
      gene_names = pp$genes, miRNA_names = pp$mirnas,
      x_center = pp$x_center, x_scale = pp$x_scale, x_median = pp$x_median,
      y_center = pp$y_center, y_scale = pp$y_scale, y_median = pp$y_median,
      score_table = score_table, selected_model = model, model_metadata = model_meta[[model]],
      selected_oof = selected_oof[[model]], split_seed = seed,
      provenance = list(source = "data/data.rda", split = "seed-specific 222/72/72 patient-level OS_status stratified split",
        calibration_role = "reliability scoring only", test_role = "final evaluation only"))
    # Correct the input-scale feature matrices after building the object.
    for (part in c("train", "cal", "test")) {
      split_name <- if (part == "cal") "calibration" else part
      ids <- partitions[[split_name]]
      xx <- dat$X[ids, pp$genes, drop = FALSE]
      for (j in seq_len(ncol(xx))) xx[!is.finite(xx[, j]), j] <- pp$x_median[j]
      canonical[[paste0("X_", part, "_native")]] <- xx
      rownames(canonical[[paste0("clinical_", part)]]) <- dat$ids[ids]
    }
    canonical$gene_names <- colnames(canonical$X_train_z)
    canonical$miRNA_names <- colnames(canonical$Y_train_z)
    can_path <- file.path(out_root, sprintf("seed%d_canonical_%s.rds", seed, model))
    saveRDS(canonical, can_path, compress = "gzip")
    canonical_paths[model] <- can_path
    log_progress(seed, model, model_meta[[model]]$config$config, "canonical_export_complete", status = "complete", elapsed_seconds = proc.time()[["elapsed"]] - t0)
  }
  # k is a PCA rank, hence not applicable to ElasticNet and Mean.
  result <- list(status = "complete", protocol_version = PROTOCOL_VERSION,
    lambda_path_policy = list(nlambda = N_LAMBDA, min_ratio_by_shape = LAMBDA_MIN_RATIO_POLICY),
    seed = seed, split = lapply(partitions, function(i) dat$ids[i]),
    split_counts = vapply(partitions, length, integer(1)), fold_ids = lapply(folds, function(i) dat$ids[tr[i]]),
    tuning = tuning, selected = aggregate, selected_model = winner, selected_oof = selected_oof, model_meta = model_meta,
    calibration_reliability = reliability, test_metrics = test_metrics, predictions = pred,
    canonical_paths = canonical_paths,
    truth = list(calibration = cap$Y_native, test = tep$Y_native),
    preprocessing = pp, k_rank_applicability = k_rank_applicability,
    k_rank_applicable = c(ElasticNet = FALSE, PCA_ElasticNet = TRUE, PCA_Ridge = TRUE, Mean = FALSE),
    elapsed_seconds = unname(proc.time()[["elapsed"]] - t0), session = capture.output(sessionInfo()))
  tmp <- paste0(target, ".tmp")
  saveRDS(result, tmp, compress = FALSE); file.rename(tmp, target)
  log_progress(seed, foldstage = "seed_summary_and_all_canonical_complete", status = "complete", elapsed_seconds = result$elapsed_seconds)
  invisible(result)
}

dat <- load_inputs(data_path)
message("Matched patients: ", length(dat$ids), "; mRNA features: ", ncol(dat$X), "; miRNA features: ", ncol(dat$Y))
if (identical(seed_arg, "sanity")) {
  p <- stratified_partition(dat$clinical$OS_status, c(train = 222L, calibration = 72L, test = 72L), 20261014L)
  stopifnot(identical(unname(vapply(p, length, integer(1))), c(222L, 72L, 72L)))
  f <- make_folds(dat$clinical$OS_status[p$train], 5L, 30261014L)
  stopifnot(length(f) == 5L, sum(lengths(f)) == 222L, !anyDuplicated(unlist(f)))
  small_x <- dat$X[p$train[1:80], 2:61, drop = FALSE]
  small_y <- dat$Y[p$train[1:80], 2:11, drop = FALSE]
  pp <- fit_preprocess(small_x[1:60, , drop = FALSE], small_y[1:60, , drop = FALSE], top_mrna = 30L)
  a <- apply_preprocess(pp, small_x[1:60, , drop = FALSE], small_y[1:60, , drop = FALSE])
  b <- apply_preprocess(pp, small_x[61:80, , drop = FALSE], small_y[61:80, , drop = FALSE])
  fit <- fit_glmnet_path(a$X, a$Y, alpha = 0.3)
  pred <- predict_native(fit, b$X, fit$lambda[5L], pp, rownames(b$X))
  stopifnot(identical(dim(pred), dim(b$Y_native)), all(is.finite(pred)))
  ms <- metric_summary(b$Y_native, pred, scale = matrix(rep(pp$y_scale, each = nrow(b$Y_native)), nrow = nrow(b$Y_native)))
  stopifnot("normalized_rmse" %in% names(ms), is.finite(ms[["normalized_rmse"]]))
  err_check <- matrix(c(2, 4, 6, 8, 10, 12), nrow = 3L, ncol = 2L)
  pred_check <- matrix(0, nrow = 3L, ncol = 2L)
  exact_vector <- sqrt(mean(c((2/2)^2, (4/2)^2, (6/2)^2, (8/4)^2, (10/4)^2, (12/4)^2)))
  got_vector <- metric_summary(err_check, pred_check, scale = c(2, 4))["normalized_rmse"]
  scale_matrix <- matrix(rep(c(2, 4), each = 3L), nrow = 3L, ncol = 2L)
  got_matrix <- metric_summary(err_check, pred_check, scale = scale_matrix)["normalized_rmse"]
  stopifnot(isTRUE(all.equal(as.numeric(got_vector), exact_vector)),
            isTRUE(all.equal(as.numeric(got_matrix), exact_vector)))
  cfs <- coef(fit, s = fit$lambda[5L]); cfm <- do.call(cbind, lapply(cfs, as.matrix))
  stopifnot(nrow(cfm) == ncol(a$X) + 1L, ncol(cfm) == ncol(a$Y))
  cat("SANITY_OK: split 222/72/72; five folds; fold-local preprocessing; multiresponse glmnet; native-scale prediction\n")
  quit(save = "no", status = 0L)
}
if (identical(seed_arg, "integration-sanity")) {
  # Full run_seed integration check on all 366 real patients but a deliberately
  # small feature subset. Its outputs are isolated from formal E1 seed files.
  sanity_out <- Sys.getenv("E1_SANITY_OUT", file.path(root, "new", "results", "_sanity", "E1"))
  out_root <- normalizePath(sanity_out, winslash = "/", mustWork = FALSE)
  if (!grepl("/new/results/_sanity(/|$)", paste0(out_root, "/"))) stop("integration-sanity output must remain under new/results/_sanity")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  progress_path <- file.path(out_root, "current_progress.csv")
  if (length(list.files(out_root, all.files = TRUE, no.. = TRUE))) stop("Sanity output directory is not empty: ", out_root)
  dat$X <- dat$X[, seq_len(min(120L, ncol(dat$X))), drop = FALSE]
  dat$Y <- dat$Y[, seq_len(min(20L, ncol(dat$Y))), drop = FALSE]
  TOP_MRNA <- 100L
  K_GRID <- c(20L, 50L)
  ALPHA_GRID <- c(0.1, 0.5)
  ans <- run_seed(20261014L, dat)
  stopifnot(ans$status == "complete", identical(ans$protocol_version, PROTOCOL_VERSION),
    all(ans$split_counts == c(train = 222L, calibration = 72L, test = 72L)),
    all(ans$tuning$nlambda_requested == N_LAMBDA),
    all(ans$tuning$lambda_min_ratio_requested_by_shape == "1e-04"))
  train_ids <- unlist(ans$split$train, use.names = FALSE)
  cal_ids <- unlist(ans$split$calibration, use.names = FALSE)
  test_ids <- unlist(ans$split$test, use.names = FALSE)
  stopifnot(!length(intersect(train_ids, cal_ids)), !length(intersect(train_ids, test_ids)), !length(intersect(cal_ids, test_ids)))
  stopifnot(setequal(unlist(lapply(ans$selected_oof, function(z) rownames(z$pred)), use.names = FALSE), rep(train_ids, times = 3L)))
  for (model in c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")) {
    can <- readRDS(ans$canonical_paths[[model]])
    stopifnot(identical(can$protocol_version, PROTOCOL_VERSION))
    stopifnot(nrow(can$X_train_z) == 222L, nrow(can$X_cal_z) == 72L, nrow(can$X_test_z) == 72L,
      nrow(can$clinical_train) == 222L, nrow(can$clinical_cal) == 72L, nrow(can$clinical_test) == 72L,
      identical(rownames(can$X_cal_z), cal_ids), identical(rownames(can$pred_cal_native), cal_ids),
      identical(as.character(can$clinical_cal$patient_id), cal_ids),
      identical(rownames(can$pred_test_native), test_ids), identical(as.character(can$clinical_test$patient_id), test_ids),
      nrow(can$score_table) == ncol(can$Y_cal_z), all(c("R2", "Pearson", "RMSE", "Rel") %in% names(can$score_table)))
    stopifnot(identical(rownames(can$selected_oof$pred), train_ids),
      length(can$selected_oof$normalized_rmse_by_miRNA$OOF_RMSE_normalized) == ncol(can$Y_train_z))
  }
  for (model in c("ElasticNet", "PCA_ElasticNet", "PCA_Ridge")) {
    stopifnot(ans$model_meta[[model]]$lambda_path$nlambda_requested == N_LAMBDA,
      ans$model_meta[[model]]$lambda_path$lambda_min_ratio_requested == 0.0001)
  }
  stopifnot(isTRUE(all.equal(ans$model_meta$PCA_Ridge$lambda_path$lambda_min_ratio_observed, 0.0001, tolerance = 1e-6)))
  # Directly assert both glmnet default lambda ranges used by the protocol.
  set.seed(40721L)
  wide_fit <- fit_glmnet_path(matrix(rnorm(12L * 20L), 12L, 20L), matrix(rnorm(12L * 3L), 12L, 3L), 0.3)
  tall_fit <- fit_glmnet_path(matrix(rnorm(20L * 12L), 20L, 12L), matrix(rnorm(20L * 3L), 20L, 3L), 0.3)
  stopifnot(length(wide_fit$lambda) <= N_LAMBDA, length(wide_fit$lambda) >= 5L,
    length(tall_fit$lambda) <= N_LAMBDA, length(tall_fit$lambda) >= 5L,
    tail(wide_fit$lambda, 1L) / wide_fit$lambda[1L] >= 0.01 * (1 - 1e-6),
    tail(tall_fit$lambda, 1L) / tall_fit$lambda[1L] >= 0.0001 * (1 - 1e-6))
  message("INTEGRATION_SANITY_OK: complete run_seed; 3 canonical model files; split IDs; OOF/cal/test provenance")
  quit(save = "no", status = 0L)
}
for (seed in seeds) {
  message("Starting seed ", seed)
  ans <- run_seed(seed, dat)
  message("Completed seed ", seed, " in ", round(ans$elapsed_seconds / 60, 1), " min; selected ", ans$selected_model)
}
