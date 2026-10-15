set.seed(123)
options(stringsAsFactors = FALSE)

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_path <- file.path(project_dir, "data", "data.rda")
out_dir <- file.path(project_dir, "output")
pathway_dir <- file.path(out_dir, "pathway_ridge_only")
dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

cache_dir <- file.path(project_dir, ".r_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
Sys.setenv(
  XDG_CACHE_HOME = cache_dir,
  R_USER_CACHE_DIR = cache_dir
)

pkgs <- c("tidyverse", "glmnet", "msigdbr", "jsonlite")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("请先安装 R 包：", pkg)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(msigdbr)
  library(jsonlite)
})

to_prediction_matrix <- function(pred, sample_names, response_names) {
  if (is.list(pred) && !is.data.frame(pred)) {
    pred <- do.call(cbind, lapply(pred, function(x) as.numeric(x)))
  }
  if (length(dim(pred)) == 3) {
    d <- dim(pred)
    if (d[1] == length(sample_names) && d[2] == length(response_names)) {
      pred <- pred[, , 1]
    } else if (d[1] == length(sample_names) && d[3] == length(response_names)) {
      pred <- pred[, 1, ]
    } else if (d[2] == length(sample_names) && d[1] == length(response_names)) {
      pred <- t(pred[, , 1])
    } else {
      pred <- matrix(as.numeric(pred), nrow = length(sample_names), ncol = length(response_names))
    }
  }
  pred <- as.matrix(pred)
  if (nrow(pred) != length(sample_names) && ncol(pred) == length(sample_names)) {
    pred <- t(pred)
  }
  rownames(pred) <- sample_names
  colnames(pred) <- response_names
  storage.mode(pred) <- "numeric"
  pred
}

build_metric_df <- function(real, pred, dataset, model_name) {
  rss <- colSums((real - pred)^2, na.rm = TRUE)
  tss <- colSums(sweep(real, 2, colMeans(real, na.rm = TRUE), "-")^2, na.rm = TRUE)
  r2 <- 1 - rss / tss
  r2[!is.finite(r2)] <- NA_real_

  pearson <- rep(NA_real_, ncol(real))
  for (j in seq_len(ncol(real))) {
    if (sd(real[, j], na.rm = TRUE) > 1e-8 && sd(pred[, j], na.rm = TRUE) > 1e-8) {
      pearson[j] <- suppressWarnings(cor(real[, j], pred[, j], use = "pairwise.complete.obs"))
    }
  }

  rmse <- sqrt(colMeans((real - pred)^2, na.rm = TRUE))

  tibble(
    dataset = dataset,
    model = model_name,
    miRNA = colnames(real),
    R2 = as.numeric(r2),
    Pearson = as.numeric(pearson),
    RMSE = as.numeric(rmse)
  )
}

build_global_df <- function(metric_df) {
  metric_df %>%
    summarise(
      model = first(model),
      mean_R2 = mean(R2, na.rm = TRUE),
      median_R2 = median(R2, na.rm = TRUE),
      mean_Pearson = mean(Pearson, na.rm = TRUE),
      median_Pearson = median(Pearson, na.rm = TRUE),
      mean_RMSE = mean(RMSE, na.rm = TRUE),
      median_RMSE = median(RMSE, na.rm = TRUE),
      n_miRNA = n()
    )
}

load(data_path)

id_col <- "patient_id"

clinical <- clinical %>%
  mutate(patient_id = as.character(.data[[id_col]])) %>%
  distinct(patient_id, .keep_all = TRUE)
mrna <- mrna %>%
  mutate(patient_id = as.character(.data[[id_col]])) %>%
  distinct(patient_id, .keep_all = TRUE)
mirna <- mirna %>%
  mutate(patient_id = as.character(.data[[id_col]])) %>%
  distinct(patient_id, .keep_all = TRUE)

common_id <- Reduce(intersect, list(clinical$patient_id, mrna$patient_id, mirna$patient_id))
mrna2 <- mrna[match(common_id, mrna$patient_id), ]
mirna2 <- mirna[match(common_id, mirna$patient_id), ]

X_expr <- mrna2[, setdiff(colnames(mrna2), id_col), drop = FALSE]
colnames(X_expr) <- make.unique(colnames(X_expr))
X_expr[] <- lapply(X_expr, function(x) suppressWarnings(as.numeric(as.character(x))))
X_raw_all <- as.matrix(X_expr)
rownames(X_raw_all) <- common_id
storage.mode(X_raw_all) <- "numeric"

Y_expr <- mirna2[, setdiff(colnames(mirna2), id_col), drop = FALSE]
colnames(Y_expr) <- make.unique(colnames(Y_expr))
Y_expr[] <- lapply(Y_expr, function(x) suppressWarnings(as.numeric(as.character(x))))
Y_raw_all <- as.matrix(Y_expr)
rownames(Y_raw_all) <- common_id
storage.mode(Y_raw_all) <- "numeric"

split_table <- readr::read_csv(file.path(out_dir, "01_train_validation_test_patient_id.csv"), show_col_types = FALSE)
train_ids <- split_table %>% filter(split == "train") %>% pull(patient_id)
val_ids <- split_table %>% filter(split == "validation") %>% pull(patient_id)
test_ids <- split_table %>% filter(split == "test") %>% pull(patient_id)

X_train_raw0 <- X_raw_all[train_ids, , drop = FALSE]
X_val_raw0 <- X_raw_all[val_ids, , drop = FALSE]
X_test_raw0 <- X_raw_all[test_ids, , drop = FALSE]
Y_train_raw0 <- Y_raw_all[train_ids, , drop = FALSE]
Y_val_raw0 <- Y_raw_all[val_ids, , drop = FALSE]
Y_test_raw0 <- Y_raw_all[test_ids, , drop = FALSE]

keep_x <- colMeans(is.finite(X_train_raw0)) >= 0.8
keep_y <- colMeans(is.finite(Y_train_raw0)) >= 0.8

X_train_raw0 <- X_train_raw0[, keep_x, drop = FALSE]
X_val_raw0 <- X_val_raw0[, keep_x, drop = FALSE]
X_test_raw0 <- X_test_raw0[, keep_x, drop = FALSE]
Y_train_raw0 <- Y_train_raw0[, keep_y, drop = FALSE]
Y_val_raw0 <- Y_val_raw0[, keep_y, drop = FALSE]
Y_test_raw0 <- Y_test_raw0[, keep_y, drop = FALSE]

x_median <- apply(X_train_raw0, 2, function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) 0 else median(x, na.rm = TRUE)
})
y_median <- apply(Y_train_raw0, 2, function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) 0 else median(x, na.rm = TRUE)
})
x_median[!is.finite(x_median)] <- 0
y_median[!is.finite(y_median)] <- 0

X_train_raw_all <- X_train_raw0
X_val_raw_all <- X_val_raw0
X_test_raw_all <- X_test_raw0
Y_train_raw <- Y_train_raw0
Y_val_raw <- Y_val_raw0
Y_test_raw <- Y_test_raw0

fill_missing <- function(mat, med) {
  for (j in seq_len(ncol(mat))) {
    idx <- !is.finite(mat[, j]) | is.na(mat[, j])
    if (any(idx)) mat[idx, j] <- med[j]
  }
  mat
}

X_train_raw_all <- fill_missing(X_train_raw_all, x_median)
X_val_raw_all <- fill_missing(X_val_raw_all, x_median)
X_test_raw_all <- fill_missing(X_test_raw_all, x_median)
Y_train_raw <- fill_missing(Y_train_raw, y_median)
Y_val_raw <- fill_missing(Y_val_raw, y_median)
Y_test_raw <- fill_missing(Y_test_raw, y_median)

x_var_all <- apply(X_train_raw_all, 2, var, na.rm = TRUE)
y_var_all <- apply(Y_train_raw, 2, var, na.rm = TRUE)
x_keep_var <- is.finite(x_var_all) & x_var_all > 1e-8
y_keep_var <- is.finite(y_var_all) & y_var_all > 1e-8

X_train_raw_all <- X_train_raw_all[, x_keep_var, drop = FALSE]
X_val_raw_all <- X_val_raw_all[, x_keep_var, drop = FALSE]
X_test_raw_all <- X_test_raw_all[, x_keep_var, drop = FALSE]
Y_train_raw <- Y_train_raw[, y_keep_var, drop = FALSE]
Y_val_raw <- Y_val_raw[, y_keep_var, drop = FALSE]
Y_test_raw <- Y_test_raw[, y_keep_var, drop = FALSE]

msig <- tryCatch(
  msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
)

gene_sets <- split(msig$gene_symbol, msig$gs_name)
gene_sets <- lapply(gene_sets, unique)

x_all_center <- colMeans(X_train_raw_all, na.rm = TRUE)
x_all_scale <- apply(X_train_raw_all, 2, sd, na.rm = TRUE)
x_all_scale[!is.finite(x_all_scale) | x_all_scale < 1e-8] <- 1

X_train_all_z <- sweep(X_train_raw_all, 2, x_all_center, "-")
X_train_all_z <- sweep(X_train_all_z, 2, x_all_scale, "/")
X_train_all_z[!is.finite(X_train_all_z)] <- 0
X_val_all_z <- sweep(X_val_raw_all, 2, x_all_center, "-")
X_val_all_z <- sweep(X_val_all_z, 2, x_all_scale, "/")
X_val_all_z[!is.finite(X_val_all_z)] <- 0
X_test_all_z <- sweep(X_test_raw_all, 2, x_all_center, "-")
X_test_all_z <- sweep(X_test_all_z, 2, x_all_scale, "/")
X_test_all_z[!is.finite(X_test_all_z)] <- 0

gene_sets_used <- lapply(gene_sets, function(gs) intersect(gs, colnames(X_train_all_z)))
gene_sets_used <- gene_sets_used[lengths(gene_sets_used) >= 5 & lengths(gene_sets_used) <= 500]
if (length(gene_sets_used) == 0) {
  stop("没有可用 pathway，请检查 mRNA 列名是否为 gene symbol")
}

build_pathway_matrix <- function(xmat, sets) {
  out <- matrix(NA_real_, nrow = nrow(xmat), ncol = length(sets), dimnames = list(rownames(xmat), names(sets)))
  for (i in seq_along(sets)) {
    out[, i] <- rowMeans(xmat[, sets[[i]], drop = FALSE], na.rm = TRUE)
  }
  out[!is.finite(out)] <- 0
  out
}

P_train_raw <- build_pathway_matrix(X_train_all_z, gene_sets_used)
P_val_raw <- build_pathway_matrix(X_val_all_z, gene_sets_used)
P_test_raw <- build_pathway_matrix(X_test_all_z, gene_sets_used)

pathway_center <- colMeans(P_train_raw, na.rm = TRUE)
pathway_scale <- apply(P_train_raw, 2, sd, na.rm = TRUE)
pathway_scale[!is.finite(pathway_scale) | pathway_scale < 1e-8] <- 1

P_train <- sweep(P_train_raw, 2, pathway_center, "-")
P_train <- sweep(P_train, 2, pathway_scale, "/")
P_train[!is.finite(P_train)] <- 0
P_val <- sweep(P_val_raw, 2, pathway_center, "-")
P_val <- sweep(P_val, 2, pathway_scale, "/")
P_val[!is.finite(P_val)] <- 0
P_test <- sweep(P_test_raw, 2, pathway_center, "-")
P_test <- sweep(P_test, 2, pathway_scale, "/")
P_test[!is.finite(P_test)] <- 0

y_center <- colMeans(Y_train_raw, na.rm = TRUE)
y_scale <- apply(Y_train_raw, 2, sd, na.rm = TRUE)
y_scale[!is.finite(y_scale) | y_scale < 1e-8] <- 1
Y_train <- sweep(Y_train_raw, 2, y_center, "-")
Y_train <- sweep(Y_train, 2, y_scale, "/")
Y_train[!is.finite(Y_train)] <- 0
Y_val <- sweep(Y_val_raw, 2, y_center, "-")
Y_val <- sweep(Y_val, 2, y_scale, "/")
Y_val[!is.finite(Y_val)] <- 0
Y_test <- sweep(Y_test_raw, 2, y_center, "-")
Y_test <- sweep(Y_test, 2, y_scale, "/")
Y_test[!is.finite(Y_test)] <- 0

fit <- cv.glmnet(
  x = P_train,
  y = Y_train,
  family = "mgaussian",
  alpha = 0,
  nfolds = 5,
  standardize = FALSE,
  type.measure = "mse"
)

pred_val <- to_prediction_matrix(predict(fit, newx = P_val, s = "lambda.min"), rownames(P_val), colnames(Y_val))
pred_test <- to_prediction_matrix(predict(fit, newx = P_test, s = "lambda.min"), rownames(P_test), colnames(Y_test))

metric_val <- build_metric_df(Y_val, pred_val, "validation", "Pathway_Ridge")
metric_test <- build_metric_df(Y_test, pred_test, "test", "Pathway_Ridge")
global_val <- build_global_df(metric_val)
global_test <- build_global_df(metric_test)

coef_list <- coef(fit, s = "lambda.min")
pathway_coef <- sapply(coef_list, function(x) as.numeric(x[-1, 1]))
pathway_coef <- as.matrix(pathway_coef)
rownames(pathway_coef) <- colnames(P_train)
colnames(pathway_coef) <- colnames(Y_train)

pathway_summary <- tibble(
  pathway = rownames(pathway_coef),
  mean_abs_coef = rowMeans(abs(pathway_coef), na.rm = TRUE),
  max_abs_coef = apply(abs(pathway_coef), 1, max, na.rm = TRUE)
) %>% arrange(desc(mean_abs_coef))

write_csv(metric_val, file.path(pathway_dir, "Pathway_Ridge_validation_metrics_by_miRNA.csv"))
write_csv(metric_test, file.path(pathway_dir, "Pathway_Ridge_test_metrics_by_miRNA.csv"))
write_csv(global_val, file.path(pathway_dir, "Pathway_Ridge_validation_global_metrics.csv"))
write_csv(global_test, file.path(pathway_dir, "Pathway_Ridge_test_global_metrics.csv"))
write_csv(pathway_summary, file.path(pathway_dir, "08_pathway_ridge_interpretation_summary.csv"))
write_csv(
  data.frame(pathway = rownames(pathway_coef), pathway_coef, check.names = FALSE),
  file.path(pathway_dir, "09_pathway_ridge_coef_matrix.csv")
)

base_val <- read_csv(file.path(out_dir, "11_validation_global_metrics.csv"), show_col_types = FALSE)
base_test <- read_csv(file.path(out_dir, "13_test_global_metrics.csv"), show_col_types = FALSE)
bind_rows(base_val, global_val) %>%
  arrange(desc(median_Pearson), desc(mean_R2), mean_RMSE) %>%
  write_csv(file.path(pathway_dir, "combined_validation_global_metrics_with_pathway.csv"))
bind_rows(base_test, global_test) %>%
  arrange(desc(median_Pearson), desc(mean_R2), mean_RMSE) %>%
  write_csv(file.path(pathway_dir, "combined_test_global_metrics_with_pathway.csv"))

summary_json <- list(
  model_name = "Pathway_Ridge",
  n_gene_sets_used = length(gene_sets_used),
  best_lambda = fit$lambda.min,
  validation_median_Pearson = global_val$median_Pearson[[1]],
  test_median_Pearson = global_test$median_Pearson[[1]],
  test_mean_R2 = global_test$mean_R2[[1]],
  test_mean_RMSE = global_test$mean_RMSE[[1]]
)

writeLines(jsonlite::toJSON(summary_json, auto_unbox = TRUE, pretty = TRUE), file.path(pathway_dir, "Pathway_Ridge_summary.json"))
cat(jsonlite::toJSON(summary_json, auto_unbox = TRUE))
