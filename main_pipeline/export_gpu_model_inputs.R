#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

data_path <- "data/data.rda"
out_dir <- "output"
gpu_dir <- file.path(out_dir, "gpu_model_inputs")
dir.create(gpu_dir, showWarnings = FALSE, recursive = TRUE)

id_col <- "patient_id"

if (!file.exists(data_path)) stop("没有找到 data/data.rda")
if (!file.exists(file.path(out_dir, "01_train_validation_test_patient_id.csv"))) {
  stop("没有找到 output/01_train_validation_test_patient_id.csv")
}
if (!file.exists(file.path(out_dir, "04_selected_high_variable_mRNA.csv"))) {
  stop("没有找到 output/04_selected_high_variable_mRNA.csv")
}

load(data_path)

if (!exists("clinical")) stop("当前环境中没有 clinical 对象")
if (!exists("mrna")) stop("当前环境中没有 mrna 对象")
if (!exists("mirna")) stop("当前环境中没有 mirna 对象")

as_unique_df <- function(df) {
  df[[id_col]] <- as.character(df[[id_col]])
  df[!duplicated(df[[id_col]]), , drop = FALSE]
}

clinical <- as_unique_df(clinical)
mrna <- as_unique_df(mrna)
mirna <- as_unique_df(mirna)

common_id <- Reduce(intersect, list(clinical[[id_col]], mrna[[id_col]], mirna[[id_col]]))

mrna2 <- mrna[match(common_id, mrna[[id_col]]), , drop = FALSE]
mirna2 <- mirna[match(common_id, mirna[[id_col]]), , drop = FALSE]

X_expr <- mrna2[, setdiff(colnames(mrna2), id_col), drop = FALSE]
colnames(X_expr) <- make.unique(colnames(X_expr))
X_expr[] <- lapply(X_expr, function(x) suppressWarnings(as.numeric(as.character(x))))
X_raw_all <- as.matrix(X_expr)
rownames(X_raw_all) <- as.character(mrna2[[id_col]])
storage.mode(X_raw_all) <- "numeric"

Y_expr <- mirna2[, setdiff(colnames(mirna2), id_col), drop = FALSE]
colnames(Y_expr) <- make.unique(colnames(Y_expr))
Y_expr[] <- lapply(Y_expr, function(x) suppressWarnings(as.numeric(as.character(x))))
Y_raw_all <- as.matrix(Y_expr)
rownames(Y_raw_all) <- as.character(mirna2[[id_col]])
storage.mode(Y_raw_all) <- "numeric"

split_table <- read.csv(file.path(out_dir, "01_train_validation_test_patient_id.csv"), check.names = FALSE)
split_table$patient_id <- as.character(split_table$patient_id)
split_table$split <- as.character(split_table$split)

split_table <- split_table[split_table$patient_id %in% rownames(X_raw_all), , drop = FALSE]

train_ids <- split_table$patient_id[split_table$split == "train"]
val_ids <- split_table$patient_id[split_table$split == "validation"]
test_ids <- split_table$patient_id[split_table$split == "test"]

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
  if (length(x) == 0) return(0)
  median(x, na.rm = TRUE)
})
y_median <- apply(Y_train_raw0, 2, function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(0)
  median(x, na.rm = TRUE)
})

x_median[!is.finite(x_median)] <- 0
y_median[!is.finite(y_median)] <- 0

impute_matrix <- function(mat, median_vec) {
  for (j in seq_len(ncol(mat))) {
    idx <- !is.finite(mat[, j]) | is.na(mat[, j])
    if (any(idx)) mat[idx, j] <- median_vec[j]
  }
  mat
}

X_train_raw_all <- impute_matrix(X_train_raw0, x_median)
X_val_raw_all <- impute_matrix(X_val_raw0, x_median)
X_test_raw_all <- impute_matrix(X_test_raw0, x_median)

Y_train_raw <- impute_matrix(Y_train_raw0, y_median)
Y_val_raw <- impute_matrix(Y_val_raw0, y_median)
Y_test_raw <- impute_matrix(Y_test_raw0, y_median)

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

top_gene_df <- read.csv(file.path(out_dir, "04_selected_high_variable_mRNA.csv"), check.names = FALSE)
top_genes <- as.character(top_gene_df$gene)
top_genes <- intersect(top_genes, colnames(X_train_raw_all))

X_train_raw <- X_train_raw_all[, top_genes, drop = FALSE]
X_val_raw <- X_val_raw_all[, top_genes, drop = FALSE]
X_test_raw <- X_test_raw_all[, top_genes, drop = FALSE]

x_center <- colMeans(X_train_raw, na.rm = TRUE)
x_scale <- apply(X_train_raw, 2, sd, na.rm = TRUE)
x_scale[!is.finite(x_scale) | x_scale < 1e-8] <- 1

y_center <- colMeans(Y_train_raw, na.rm = TRUE)
y_scale <- apply(Y_train_raw, 2, sd, na.rm = TRUE)
y_scale[!is.finite(y_scale) | y_scale < 1e-8] <- 1

scale_matrix <- function(mat, center_vec, scale_vec) {
  out <- sweep(mat, 2, center_vec, "-")
  out <- sweep(out, 2, scale_vec, "/")
  out[!is.finite(out)] <- 0
  out
}

X_train <- scale_matrix(X_train_raw, x_center, x_scale)
X_val <- scale_matrix(X_val_raw, x_center, x_scale)
X_test <- scale_matrix(X_test_raw, x_center, x_scale)

Y_train <- scale_matrix(Y_train_raw, y_center, y_scale)
Y_val <- scale_matrix(Y_val_raw, y_center, y_scale)
Y_test <- scale_matrix(Y_test_raw, y_center, y_scale)

write_with_id <- function(mat, path) {
  df <- data.frame(patient_id = rownames(mat), mat, check.names = FALSE)
  write.csv(df, path, row.names = FALSE)
}

write_with_id(X_train, file.path(gpu_dir, "train_x_z.csv"))
write_with_id(X_val, file.path(gpu_dir, "validation_x_z.csv"))
write_with_id(X_test, file.path(gpu_dir, "test_x_z.csv"))

write_with_id(Y_train, file.path(gpu_dir, "train_y_z.csv"))
write_with_id(Y_val, file.path(gpu_dir, "validation_y_z.csv"))
write_with_id(Y_test, file.path(gpu_dir, "test_y_z.csv"))

write.csv(
  data.frame(
    feature = names(x_center),
    x_center = as.numeric(x_center),
    x_scale = as.numeric(x_scale),
    check.names = FALSE
  ),
  file.path(gpu_dir, "x_scaling_stats.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    miRNA = names(y_center),
    y_center = as.numeric(y_center),
    y_scale = as.numeric(y_scale),
    check.names = FALSE
  ),
  file.path(gpu_dir, "y_scaling_stats.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    split = c("train", "validation", "test"),
    n = c(nrow(X_train), nrow(X_val), nrow(X_test)),
    p = c(ncol(X_train), ncol(X_val), ncol(X_test)),
    q = c(ncol(Y_train), ncol(Y_val), ncol(Y_test))
  ),
  file.path(gpu_dir, "gpu_input_summary.csv"),
  row.names = FALSE
)

message("GPU 模型输入已导出到：", gpu_dir)
