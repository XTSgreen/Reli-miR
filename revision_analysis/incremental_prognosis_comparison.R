# Three-model prognostic comparison on the frozen E0 fixed split.
# Tests whether virtual miRNA (fixed PCA_ElasticNet predictions) adds
# incremental prognostic value beyond mRNA input (pre-registered:
# new/results/E5/PREREGISTRATION_threemodel_comparison.md).
suppressWarnings(suppressMessages({
  library(survival); library(glmnet)
}))
B <- 2000L
BOOT_SEED <- 20261019L
ALPHA <- 0.5
NFOLDS <- 5L

source("code/run_external_generalization.R", local = TRUE)  # to_prediction_matrix helper

x <- readRDS("new/results/E0/canonical_input.rds")
obj <- x$metadata_lightweight_object
if (is.null(obj)) obj <- readRDS(x$lightweight_object_path)
model <- readRDS(x$metadata_model_path)
stopifnot(!is.null(model$pca), !is.null(model$fit))

predict_virtual_z <- function(X_z) {
  pc <- predict(model$pca, newdata = X_z[, rownames(model$pca$rotation), drop = FALSE])
  raw <- predict(model$fit, newx = pc[, seq_len(model$best_k), drop = FALSE], s = "lambda.min")
  to_prediction_matrix(raw, rownames(X_z), colnames(x$Y_test_z))
}

clin <- function(d) {
  time <- suppressWarnings(as.numeric(as.character(d$OS_time)))
  status <- suppressWarnings(as.integer(as.character(d$OS_status)))
  data.frame(patient_id = d$patient_id, time = time, status = status)
}
tr_cl <- clin(x$clinical_train); te_cl <- clin(x$clinical_test)
stopifnot(identical(tr_cl$patient_id, rownames(x$X_train_z)),
          identical(te_cl$patient_id, rownames(x$X_test_z)))

M_tr <- as.matrix(x$X_train_z); M_te <- as.matrix(x$X_test_z)
V_tr <- predict_virtual_z(M_tr); V_te <- predict_virtual_z(M_te)
colnames(V_tr) <- paste0("virt__", colnames(V_tr)); colnames(V_te) <- colnames(V_tr)
R_tr <- as.matrix(x$Y_train_z); R_te <- as.matrix(x$Y_test_z)
colnames(R_tr) <- paste0("real__", colnames(R_tr)); colnames(R_te) <- colnames(R_tr)

# Pre-registered amendment: drop training rows with OS_time <= 0 (2 rows, 0 events)
# before constructing any feature set; glmnet Cox forbids non-positive times.
keep_tr <- is.finite(tr_cl$time) & tr_cl$time > 0
M_tr <- M_tr[keep_tr, , drop = FALSE]
V_tr <- V_tr[keep_tr, , drop = FALSE]
R_tr <- R_tr[keep_tr, , drop = FALSE]
y_tr <- Surv(tr_cl$time, tr_cl$status)[keep_tr, ]
cat("training rows kept:", sum(keep_tr), "of", length(keep_tr),
    "; events:", sum(tr_cl$status[keep_tr] == 1L), "\n")

sets <- list(
  M_mRNA = list(tr = M_tr, te = M_te),
  V_virtual = list(tr = V_tr, te = V_te),
  C_combined = list(tr = cbind(M_tr, V_tr), te = cbind(M_te, V_te)),
  R_real_diag = list(tr = R_tr, te = R_te))

fit_cox <- function(X, y, seed) {
  # Fixed per-model CV fold seed for reproducibility (pre-registration amendment 9).
  set.seed(seed)
  cvfit <- cv.glmnet(X, y, family = "cox", alpha = ALPHA, nfolds = NFOLDS,
                     standardize = FALSE, cox.ties = "efron")
  coef_vec <- as.matrix(coef(cvfit, s = "lambda.min"))
  lp <- drop(X %*% coef_vec)
  list(cvfit = cvfit, coef = coef_vec, lp = lp,
       n_active = sum(coef_vec != 0))
}

fits <- list()
for (i in seq_along(sets)) {
  nm <- names(sets)[i]
  fits[[nm]] <- fit_cox(sets[[nm]]$tr, y_tr, seed = 20261019L + i)
  cat(nm, "active=", fits[[nm]]$n_active, " lambda.min=", fits[[nm]]$cvfit$lambda.min, "\n")
}

cindex <- function(lp, time, status) {
  ok <- is.finite(lp) & is.finite(time) & is.finite(status)
  if (sum(ok) < 5L || sum(status[ok] == 1L) < 2L) return(NA_real_)
  cc <- survival::concordance(Surv(time[ok], status[ok]) ~ lp[ok])
  # The formula interface of survival::concordance returns 1 - Harrell C for
  # Cox-type risk scores (verified against coxph summary concordance);
  # convert to the standard orientation (higher lp = higher hazard).
  1 - as.numeric(cc["concordance"])
}

te_lp <- sapply(names(sets), function(nm) drop(sets[[nm]]$te %*% fits[[nm]]$coef))
te_c <- sapply(names(sets), function(nm) cindex(te_lp[, nm], te_cl$time, te_cl$status))
# Direction sanity check: in-sample training C should exceed 0.5 if the fit
# orientation and the lp sign are correct.
tr_c_insample <- sapply(names(sets), function(nm) cindex(fits[[nm]]$lp,
  as.numeric(y_tr[, "time"]), as.numeric(y_tr[, "status"])))
cat("train in-sample C:\n"); print(round(tr_c_insample, 4))
cat("test C-index:\n"); print(round(te_c, 4))

set.seed(BOOT_SEED)
n_te <- nrow(M_te)
contrast_list <- list("V_minus_M" = c("V_virtual", "M_mRNA"),
                      "C_minus_M" = c("C_combined", "M_mRNA"),
                      "C_minus_V" = c("C_combined", "V_virtual"),
                      "R_minus_M" = c("R_real_diag", "M_mRNA"))
bs <- vector("list", B)
for (b in seq_len(B)) {
  ii <- sample.int(n_te, n_te, replace = TRUE)
  out_b <- sapply(names(sets), function(nm) cindex(te_lp[ii, nm], te_cl$time[ii], te_cl$status[ii]))
  for (nm in names(contrast_list)) {
    pr <- contrast_list[[nm]]
    out_b[[nm]] <- out_b[[pr[1]]] - out_b[[pr[2]]]
  }
  bs[[b]] <- out_b
}
bs_df <- as.data.frame(do.call(rbind, bs), stringsAsFactors = FALSE)
# point estimates for differences (same patients, no resampling)
diff_points <- sapply(contrast_list, function(pr)
  cindex(te_lp[, pr[1]], te_cl$time, te_cl$status) -
    cindex(te_lp[, pr[2]], te_cl$time, te_cl$status))
summ <- do.call(rbind, lapply(names(bs_df), function(nm) {
  v <- bs_df[[nm]]
  pt <- if (nm %in% names(sets)) unname(te_c[[nm]]) else unname(diff_points[[nm]])
  data.frame(metric = nm, point_estimate = pt,
             boot_median = median(v, na.rm = TRUE),
             CI_low = quantile(v, .025, na.rm = TRUE),
             CI_high = quantile(v, .975, na.rm = TRUE),
             B_valid = sum(is.finite(v)), B_requested = B)
}))
out_dir <- "new/results/E5"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(summ, file.path(out_dir, "E5_threemodel_survival_comparison.csv"), row.names = FALSE)
meta <- data.frame(protocol = "fixed_split_E0", alpha = ALPHA, nfolds = NFOLDS,
  lambda_min = sapply(fits, function(f) f$cvfit$lambda.min),
  n_active = sapply(fits, function(f) f$n_active),
  train_insample_C = tr_c_insample[names(fits)],
  n_train = nrow(M_tr), n_test = n_te,
  train_events = sum(tr_cl$status[keep_tr] == 1L), test_events = sum(te_cl$status == 1L),
  B = B, boot_seed = BOOT_SEED)
write.csv(meta, file.path(out_dir, "E5_threemodel_run_metadata.csv"), row.names = TRUE)
cat("C-index differences:\n"); print(round(diff_points, 4))
cat("summary written\n"); print(summ)
