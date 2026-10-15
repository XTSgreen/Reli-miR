#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(glmnet)
})

# Run from the repository root; override with --project-dir if needed.
project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "output", "gpu_model_inputs")
out_dir <- file.path(project_dir, "output", "pca_elasticnet_only")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_matrix_csv <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  ids <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  mat
}

to_prediction_matrix <- function(pred, sample_names, response_names) {
  if (is.list(pred) && !is.data.frame(pred)) {
    pred <- do.call(cbind, lapply(pred, function(x) as.numeric(x)))
  }

  if (length(dim(pred)) == 3) {
    d <- dim(pred)

    if (d[1] == length(sample_names) && d[2] == length(response_names)) {
      pred <- pred[, , 1, drop = FALSE]
      pred <- pred[, , 1]
    } else if (d[1] == length(sample_names) && d[3] == length(response_names)) {
      pred <- pred[, 1, , drop = FALSE]
      pred <- pred[, 1, ]
    } else if (d[2] == length(sample_names) && d[1] == length(response_names)) {
      pred <- pred[, , 1]
      pred <- t(pred)
    } else {
      pred <- matrix(as.numeric(pred), nrow = length(sample_names), ncol = length(response_names))
    }
  }

  pred <- as.matrix(pred)

  if (nrow(pred) != length(sample_names) && ncol(pred) == length(sample_names)) {
    pred <- t(pred)
  }

  if (ncol(pred) != length(response_names)) {
    pred <- matrix(as.numeric(pred), nrow = length(sample_names), ncol = length(response_names))
  }

  rownames(pred) <- sample_names
  colnames(pred) <- response_names
  storage.mode(pred) <- "numeric"
  pred
}

build_metric_df <- function(dataset, model_name, real, pred) {
  rss <- colSums((real - pred)^2, na.rm = TRUE)
  tss <- colSums(sweep(real, 2, colMeans(real, na.rm = TRUE), "-")^2, na.rm = TRUE)

  r2 <- 1 - rss / tss
  r2[!is.finite(r2)] <- NA_real_

  pearson <- rep(NA_real_, ncol(real))
  for (j in seq_len(ncol(real))) {
    if (sd(real[, j], na.rm = TRUE) > 1e-8 && sd(pred[, j], na.rm = TRUE) > 1e-8) {
      pearson[j] <- suppressWarnings(cor(real[, j], pred[, j], use = "complete.obs"))
    }
  }

  rmse <- sqrt(colMeans((real - pred)^2, na.rm = TRUE))

  data.frame(
    dataset = dataset,
    model = model_name,
    miRNA = colnames(real),
    R2 = r2,
    Pearson = pearson,
    RMSE = rmse,
    check.names = FALSE
  )
}

summarise_metric_df <- function(metric_df) {
  data.frame(
    model = unique(metric_df$model),
    mean_R2 = mean(metric_df$R2, na.rm = TRUE),
    median_R2 = median(metric_df$R2, na.rm = TRUE),
    mean_Pearson = mean(metric_df$Pearson, na.rm = TRUE),
    median_Pearson = median(metric_df$Pearson, na.rm = TRUE),
    mean_RMSE = mean(metric_df$RMSE, na.rm = TRUE),
    median_RMSE = median(metric_df$RMSE, na.rm = TRUE),
    n_miRNA = nrow(metric_df),
    check.names = FALSE
  )
}

X_train <- read_matrix_csv(file.path(input_dir, "train_x_z.csv"))
X_val <- read_matrix_csv(file.path(input_dir, "validation_x_z.csv"))
X_test <- read_matrix_csv(file.path(input_dir, "test_x_z.csv"))

Y_train <- read_matrix_csv(file.path(input_dir, "train_y_z.csv"))
Y_val <- read_matrix_csv(file.path(input_dir, "validation_y_z.csv"))
Y_test <- read_matrix_csv(file.path(input_dir, "test_y_z.csv"))

miRNA_names <- colnames(Y_train)
nfolds <- 5
k_grid <- c(20, 50, 100, 150, 200)
alpha_grid <- c(0.1, 0.3, 0.5, 0.7, 0.9)

max_k <- min(nrow(X_train) - 2, ncol(X_train), max(k_grid))
k_grid <- k_grid[k_grid <= max_k]

if (length(k_grid) == 0) {
  stop("PCA-ElasticNet 可用 PC 数不足")
}

message("Training PCA-ElasticNet only runner")

pca_fit <- prcomp(X_train, center = FALSE, scale. = FALSE)

tuning_rows <- list()
best_fit <- NULL
best_k <- NA_real_
best_alpha <- NA_real_
best_rmse <- Inf

for (k in k_grid) {
  message("PCA-ElasticNet k = ", k)
  Z_train <- pca_fit$x[, 1:k, drop = FALSE]
  Z_val <- predict(pca_fit, newdata = X_val)[, 1:k, drop = FALSE]

  for (alpha in alpha_grid) {
    message("  alpha = ", alpha)

    fit_ka <- cv.glmnet(
      x = Z_train,
      y = Y_train,
      family = "mgaussian",
      alpha = alpha,
      nfolds = nfolds,
      standardize = FALSE,
      type.measure = "mse"
    )

    pred_val <- predict(fit_ka, newx = Z_val, s = "lambda.min")
    pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)
    val_rmse <- sqrt(mean((Y_val - pred_val)^2, na.rm = TRUE))

    tuning_rows[[length(tuning_rows) + 1]] <- data.frame(
      k = k,
      alpha = alpha,
      validation_RMSE = val_rmse
    )

    if (is.finite(val_rmse) && val_rmse < best_rmse) {
      best_rmse <- val_rmse
      best_fit <- fit_ka
      best_k <- k
      best_alpha <- alpha
    }
  }
}

if (is.null(best_fit) || !is.finite(best_rmse)) {
  stop("PCA-ElasticNet 没有找到有效超参数组合")
}

Z_val <- predict(pca_fit, newdata = X_val)[, 1:best_k, drop = FALSE]
Z_test <- predict(pca_fit, newdata = X_test)[, 1:best_k, drop = FALSE]

pred_val <- predict(best_fit, newx = Z_val, s = "lambda.min")
pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)

pred_test <- predict(best_fit, newx = Z_test, s = "lambda.min")
pred_test <- to_prediction_matrix(pred_test, rownames(X_test), miRNA_names)

metric_val <- build_metric_df("validation", "PCA_ElasticNet", Y_val, pred_val)
metric_test <- build_metric_df("test", "PCA_ElasticNet", Y_test, pred_test)

global_val <- summarise_metric_df(metric_val)
global_test <- summarise_metric_df(metric_test)

write.csv(do.call(rbind, tuning_rows), file.path(out_dir, "07b_pca_elasticnet_selection.csv"), row.names = FALSE)
write.csv(metric_val, file.path(out_dir, "PCA_ElasticNet_validation_metrics_by_miRNA.csv"), row.names = FALSE)
write.csv(metric_test, file.path(out_dir, "PCA_ElasticNet_test_metrics_by_miRNA.csv"), row.names = FALSE)
write.csv(global_val, file.path(out_dir, "PCA_ElasticNet_validation_global_metrics.csv"), row.names = FALSE)
write.csv(global_test, file.path(out_dir, "PCA_ElasticNet_test_global_metrics.csv"), row.names = FALSE)

write.csv(
  data.frame(patient_id = rownames(pred_val), pred_val, check.names = FALSE),
  file.path(out_dir, "PCA_ElasticNet_prediction_validation_z.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(patient_id = rownames(pred_test), pred_test, check.names = FALSE),
  file.path(out_dir, "PCA_ElasticNet_prediction_test_z.csv"),
  row.names = FALSE
)

base_val <- read.csv(file.path(project_dir, "output", "11_validation_global_metrics.csv"), check.names = FALSE)
base_test <- read.csv(file.path(project_dir, "output", "13_test_global_metrics.csv"), check.names = FALSE)

combined_val <- rbind(base_val, global_val)
combined_val <- combined_val[order(-combined_val$median_Pearson, -combined_val$mean_R2, combined_val$mean_RMSE), , drop = FALSE]

combined_test <- rbind(base_test, global_test)
combined_test <- combined_test[order(-combined_test$median_Pearson, -combined_test$mean_R2, combined_test$mean_RMSE), , drop = FALSE]

write.csv(combined_val, file.path(out_dir, "combined_validation_global_metrics_with_pca_elasticnet.csv"), row.names = FALSE)
write.csv(combined_test, file.path(out_dir, "combined_test_global_metrics_with_pca_elasticnet.csv"), row.names = FALSE)

summary_df <- data.frame(
  model = "PCA_ElasticNet",
  best_k = best_k,
  best_alpha = best_alpha,
  best_validation_rmse = best_rmse,
  validation_median_Pearson = global_val$median_Pearson,
  test_median_Pearson = global_test$median_Pearson,
  check.names = FALSE
)
write.csv(summary_df, file.path(out_dir, "PCA_ElasticNet_summary.csv"), row.names = FALSE)

message("PCA_ElasticNet only runner finished.")
