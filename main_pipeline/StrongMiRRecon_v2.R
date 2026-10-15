# ============================================================
# Reli-miR final script
# 顺序脚本版：mRNA-to-miRNA virtual reconstruction framework
#
# 适用场景：
#   已有 paired mRNA-miRNA 数据集：
#     clinical: patient_id + 临床信息
#     mrna    : patient_id + mRNA expression matrix
#     mirna   : patient_id + miRNA expression matrix
#
# 核心改进：
#   1. train / validation / test 三划分
#   2. validation 选择最佳模型，test 只做最终评价
#   3. 增加 Mean baseline 和 Shuffled baseline
#   4. 增加 miRNA-level reliability
#   5. 增加 sample-level reliability
#   6. 私有 mRNA 预测时增加 gene coverage 检查
#   7. 保存模型、预测矩阵、样本划分、sessionInfo 和主要结果
#   8. 按 miRNA-level reliability 阈值分别绘制 downstream preservation 图
#
# 建议文件路径：
#   data/data.rda
#   data.rda 中包含 clinical、mrna、mirna 三个对象
# ============================================================


# 0. 基础参数 ---------------------------------------------------------------

set.seed(123)
options(stringsAsFactors = FALSE)

# 将 R 缓存固定到项目目录下，避免 Windows 用户目录缓存被沙箱拦截。
cache_dir <- Sys.getenv("STRONGMIR_R_CACHE", unset = file.path(".", ".r_cache"))
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cache_dir <- normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
Sys.setenv(
  XDG_CACHE_HOME = cache_dir,
  R_USER_CACHE_DIR = cache_dir
)

data_path <- "data/data.rda"

out_dir <- "output"
fig_dir <- file.path(out_dir, "figures")
model_dir <- file.path(out_dir, "models")
pred_dir  <- file.path(out_dir, "predictions")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pred_dir, showWarnings = FALSE, recursive = TRUE)

id_col <- "patient_id"

top_mrna <- 3000

train_prop <- 0.60
val_prop   <- 0.20
test_prop  <- 0.20

nfolds <- 5

run_mlp <- tolower(Sys.getenv("STRONGMIR_RUN_MLP", "true")) == "true"
run_pathway_ridge <- tolower(Sys.getenv("STRONGMIR_RUN_PATHWAY_RIDGE", "true")) == "true"
run_shuffle_baseline <- tolower(Sys.getenv("STRONGMIR_RUN_SHUFFLE_BASELINE", "true")) == "true"
run_pca_elasticnet <- tolower(Sys.getenv("STRONGMIR_RUN_PCA_ELASTICNET", "true")) == "true"

# 这些是真正参与最佳模型选择的模型
candidate_models <- c(
  "Ridge",
  "ElasticNet",
  "PLS",
  "PCA_Ridge",
  "PCA_ElasticNet",
  "MLP",
  "Pathway_Ridge"
)


# 1. 加载 R 包 --------------------------------------------------------------

pkgs <- c(
  "tidyverse",
  "glmnet",
  "pls",
  "survival",
  "nnet"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("请先安装 R 包：", pkg)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(pls)
  library(survival)
  library(nnet)
})


# 通用预测矩阵整理函数 -------------------------------------------------------
# 解决 glmnet::predict(family = "mgaussian") 有时返回 3D array / list，
# 导致 colnames 赋值时报 “dimnames 长度不匹配” 的问题。
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
      pred <- matrix(
        as.numeric(pred),
        nrow = length(sample_names),
        ncol = length(response_names)
      )
    }
  }

  pred <- as.matrix(pred)

  if (nrow(pred) != length(sample_names) && ncol(pred) == length(sample_names)) {
    pred <- t(pred)
  }

  if (ncol(pred) != length(response_names)) {
    if (length(pred) == length(sample_names) * length(response_names)) {
      pred <- matrix(
        as.numeric(pred),
        nrow = length(sample_names),
        ncol = length(response_names)
      )
    } else {
      stop(
        "预测矩阵维度不匹配：得到 ", nrow(pred), " x ", ncol(pred),
        "，但期望 ", length(sample_names), " x ", length(response_names)
      )
    }
  }

  rownames(pred) <- sample_names
  colnames(pred) <- response_names
  storage.mode(pred) <- "numeric"

  pred
}


# 2. 加载数据 ---------------------------------------------------------------

if (file.exists(data_path)) {
  load(data_path)
} else {
  message("没有找到 data/data.rda。")
  message("如果你已经在当前 R 环境中加载了 clinical、mrna、mirna，可以继续运行。")
}

if (!exists("clinical")) stop("当前环境中没有 clinical 对象")
if (!exists("mrna")) stop("当前环境中没有 mrna 对象")
if (!exists("mirna")) stop("当前环境中没有 mirna 对象")

if (!id_col %in% colnames(clinical)) stop("clinical 中没有 patient_id 列")
if (!id_col %in% colnames(mrna)) stop("mrna 中没有 patient_id 列")
if (!id_col %in% colnames(mirna)) stop("mirna 中没有 patient_id 列")


# 3. 数据对齐 ---------------------------------------------------------------

clinical <- clinical %>%
  mutate(patient_id = as.character(.data[[id_col]])) %>%
  distinct(patient_id, .keep_all = TRUE)

mrna <- mrna %>%
  mutate(patient_id = as.character(.data[[id_col]])) %>%
  distinct(patient_id, .keep_all = TRUE)

mirna <- mirna %>%
  mutate(patient_id = as.character(.data[[id_col]])) %>%
  distinct(patient_id, .keep_all = TRUE)

common_id <- Reduce(
  intersect,
  list(
    clinical$patient_id,
    mrna$patient_id,
    mirna$patient_id
  )
)

message("配对样本数量：", length(common_id))

if (length(common_id) < 50) {
  warning("配对样本数量较少，模型性能和 downstream 分析可能不稳定。")
}

clinical2 <- clinical[match(common_id, clinical$patient_id), ]
mrna2 <- mrna[match(common_id, mrna$patient_id), ]
mirna2 <- mirna[match(common_id, mirna$patient_id), ]

# mRNA 转换成 numeric matrix
X_ids <- as.character(mrna2[[id_col]])

X_expr <- mrna2[, setdiff(colnames(mrna2), id_col), drop = FALSE]
colnames(X_expr) <- make.unique(colnames(X_expr))

X_expr[] <- lapply(X_expr, function(x) {
  suppressWarnings(as.numeric(as.character(x)))
})

X_raw_all <- as.matrix(X_expr)
rownames(X_raw_all) <- X_ids
storage.mode(X_raw_all) <- "numeric"

# miRNA 转换成 numeric matrix
Y_ids <- as.character(mirna2[[id_col]])

Y_expr <- mirna2[, setdiff(colnames(mirna2), id_col), drop = FALSE]
colnames(Y_expr) <- make.unique(colnames(Y_expr))

Y_expr[] <- lapply(Y_expr, function(x) {
  suppressWarnings(as.numeric(as.character(x)))
})

Y_raw_all <- as.matrix(Y_expr)
rownames(Y_raw_all) <- Y_ids
storage.mode(Y_raw_all) <- "numeric"

# 确保顺序完全一致
X_raw_all <- X_raw_all[common_id, , drop = FALSE]
Y_raw_all <- Y_raw_all[common_id, , drop = FALSE]


# 4. train / validation / test 三划分 ---------------------------------------
# 说明：
#   如果存在 OS_status，则按 OS_status 分层划分
#   否则普通随机划分
#   最佳模型只用 validation 选择
#   test 只做最终评价

set.seed(123)

n <- nrow(clinical2)

train_idx <- integer(0)
val_idx <- integer(0)
test_idx <- integer(0)

if ("OS_status" %in% colnames(clinical2)) {

  split_group <- as.character(clinical2$OS_status)
  split_group[is.na(split_group)] <- "NA"

  for (g in unique(split_group)) {

    idx_g <- which(split_group == g)
    idx_g <- sample(idx_g)

    n_g <- length(idx_g)

    n_test_g <- floor(n_g * test_prop)
    n_val_g  <- floor(n_g * val_prop)

    if (n_g >= 10) {
      n_test_g <- max(1, n_test_g)
      n_val_g  <- max(1, n_val_g)
    }

    test_part <- integer(0)
    val_part <- integer(0)

    if (n_test_g > 0) {
      test_part <- idx_g[seq_len(n_test_g)]
    }

    if (n_val_g > 0) {
      val_start <- n_test_g + 1
      val_end <- n_test_g + n_val_g
      val_part <- idx_g[val_start:val_end]
    }

    train_part <- setdiff(idx_g, c(test_part, val_part))

    test_idx <- c(test_idx, test_part)
    val_idx <- c(val_idx, val_part)
    train_idx <- c(train_idx, train_part)
  }

} else {

  all_idx <- sample(seq_len(n))

  n_test <- floor(n * test_prop)
  n_val  <- floor(n * val_prop)

  test_idx <- all_idx[seq_len(n_test)]
  val_idx <- all_idx[(n_test + 1):(n_test + n_val)]
  train_idx <- setdiff(all_idx, c(test_idx, val_idx))
}

train_idx <- sort(unique(train_idx))
val_idx <- sort(unique(val_idx))
test_idx <- sort(unique(test_idx))

clinical_train <- clinical2[train_idx, , drop = FALSE]
clinical_val   <- clinical2[val_idx, , drop = FALSE]
clinical_test  <- clinical2[test_idx, , drop = FALSE]

X_train_raw0 <- X_raw_all[train_idx, , drop = FALSE]
X_val_raw0   <- X_raw_all[val_idx, , drop = FALSE]
X_test_raw0  <- X_raw_all[test_idx, , drop = FALSE]

Y_train_raw0 <- Y_raw_all[train_idx, , drop = FALSE]
Y_val_raw0   <- Y_raw_all[val_idx, , drop = FALSE]
Y_test_raw0  <- Y_raw_all[test_idx, , drop = FALSE]

split_table <- tibble(
  patient_id = c(
    rownames(X_train_raw0),
    rownames(X_val_raw0),
    rownames(X_test_raw0)
  ),
  split = c(
    rep("train", nrow(X_train_raw0)),
    rep("validation", nrow(X_val_raw0)),
    rep("test", nrow(X_test_raw0))
  )
)

write_csv(
  split_table,
  file.path(out_dir, "01_train_validation_test_patient_id.csv")
)

write_csv(
  tibble(
    split = c("train", "validation", "test"),
    n = c(nrow(X_train_raw0), nrow(X_val_raw0), nrow(X_test_raw0))
  ),
  file.path(out_dir, "02_train_validation_test_summary.csv")
)


# 5. 预处理 -----------------------------------------------------------------
# 原则：
#   所有筛选、填补、标准化参数都只从 train 中学习
#   validation 和 test 只应用 train 得到的参数

# 5.1 过滤缺失比例过高的 mRNA 和 miRNA
keep_x <- colMeans(is.finite(X_train_raw0)) >= 0.8
keep_y <- colMeans(is.finite(Y_train_raw0)) >= 0.8

X_train_raw0 <- X_train_raw0[, keep_x, drop = FALSE]
X_val_raw0   <- X_val_raw0[, keep_x, drop = FALSE]
X_test_raw0  <- X_test_raw0[, keep_x, drop = FALSE]

Y_train_raw0 <- Y_train_raw0[, keep_y, drop = FALSE]
Y_val_raw0   <- Y_val_raw0[, keep_y, drop = FALSE]
Y_test_raw0  <- Y_test_raw0[, keep_y, drop = FALSE]

# 5.2 用训练集 median 填补缺失值
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

X_train_raw_all <- X_train_raw0
X_val_raw_all   <- X_val_raw0
X_test_raw_all  <- X_test_raw0

Y_train_raw <- Y_train_raw0
Y_val_raw   <- Y_val_raw0
Y_test_raw  <- Y_test_raw0

for (j in seq_len(ncol(X_train_raw_all))) {

  idx <- !is.finite(X_train_raw_all[, j]) | is.na(X_train_raw_all[, j])
  if (any(idx)) X_train_raw_all[idx, j] <- x_median[j]

  idx <- !is.finite(X_val_raw_all[, j]) | is.na(X_val_raw_all[, j])
  if (any(idx)) X_val_raw_all[idx, j] <- x_median[j]

  idx <- !is.finite(X_test_raw_all[, j]) | is.na(X_test_raw_all[, j])
  if (any(idx)) X_test_raw_all[idx, j] <- x_median[j]
}

for (j in seq_len(ncol(Y_train_raw))) {

  idx <- !is.finite(Y_train_raw[, j]) | is.na(Y_train_raw[, j])
  if (any(idx)) Y_train_raw[idx, j] <- y_median[j]

  idx <- !is.finite(Y_val_raw[, j]) | is.na(Y_val_raw[, j])
  if (any(idx)) Y_val_raw[idx, j] <- y_median[j]

  idx <- !is.finite(Y_test_raw[, j]) | is.na(Y_test_raw[, j])
  if (any(idx)) Y_test_raw[idx, j] <- y_median[j]
}

# 5.3 去除训练集中低方差 mRNA 和 miRNA
x_var_all <- apply(X_train_raw_all, 2, var, na.rm = TRUE)
y_var_all <- apply(Y_train_raw, 2, var, na.rm = TRUE)

x_keep_var <- is.finite(x_var_all) & x_var_all > 1e-8
y_keep_var <- is.finite(y_var_all) & y_var_all > 1e-8

X_train_raw_all <- X_train_raw_all[, x_keep_var, drop = FALSE]
X_val_raw_all   <- X_val_raw_all[, x_keep_var, drop = FALSE]
X_test_raw_all  <- X_test_raw_all[, x_keep_var, drop = FALSE]

Y_train_raw <- Y_train_raw[, y_keep_var, drop = FALSE]
Y_val_raw   <- Y_val_raw[, y_keep_var, drop = FALSE]
Y_test_raw  <- Y_test_raw[, y_keep_var, drop = FALSE]

x_median <- x_median[colnames(X_train_raw_all)]
y_median <- y_median[colnames(Y_train_raw)]

# 5.4 高变异 mRNA 筛选
x_var <- apply(X_train_raw_all, 2, var, na.rm = TRUE)

top_genes <- names(sort(x_var, decreasing = TRUE))[1:min(top_mrna, length(x_var))]

X_train_raw <- X_train_raw_all[, top_genes, drop = FALSE]
X_val_raw   <- X_val_raw_all[, top_genes, drop = FALSE]
X_test_raw  <- X_test_raw_all[, top_genes, drop = FALSE]

# 5.5 Z-score 标准化
x_center <- colMeans(X_train_raw, na.rm = TRUE)
x_scale <- apply(X_train_raw, 2, sd, na.rm = TRUE)
x_scale[!is.finite(x_scale) | x_scale < 1e-8] <- 1

y_center <- colMeans(Y_train_raw, na.rm = TRUE)
y_scale <- apply(Y_train_raw, 2, sd, na.rm = TRUE)
y_scale[!is.finite(y_scale) | y_scale < 1e-8] <- 1

X_train <- sweep(X_train_raw, 2, x_center, "-")
X_train <- sweep(X_train, 2, x_scale, "/")
X_train[!is.finite(X_train)] <- 0

X_val <- sweep(X_val_raw, 2, x_center, "-")
X_val <- sweep(X_val, 2, x_scale, "/")
X_val[!is.finite(X_val)] <- 0

X_test <- sweep(X_test_raw, 2, x_center, "-")
X_test <- sweep(X_test, 2, x_scale, "/")
X_test[!is.finite(X_test)] <- 0

Y_train <- sweep(Y_train_raw, 2, y_center, "-")
Y_train <- sweep(Y_train, 2, y_scale, "/")
Y_train[!is.finite(Y_train)] <- 0

Y_val <- sweep(Y_val_raw, 2, y_center, "-")
Y_val <- sweep(Y_val, 2, y_scale, "/")
Y_val[!is.finite(Y_val)] <- 0

Y_test <- sweep(Y_test_raw, 2, y_center, "-")
Y_test <- sweep(Y_test, 2, y_scale, "/")
Y_test[!is.finite(Y_test)] <- 0

miRNA_names <- colnames(Y_train)

write_csv(
  tibble(
    item = c(
      "n_all_paired_samples",
      "n_train",
      "n_validation",
      "n_test",
      "n_mRNA_after_filter",
      "n_selected_mRNA",
      "n_miRNA_after_filter"
    ),
    value = c(
      length(common_id),
      nrow(X_train),
      nrow(X_val),
      nrow(X_test),
      ncol(X_train_raw_all),
      ncol(X_train),
      ncol(Y_train)
    )
  ),
  file.path(out_dir, "03_preprocess_summary.csv")
)

write_csv(
  tibble(
    gene = top_genes,
    variance = x_var[top_genes]
  ),
  file.path(out_dir, "04_selected_high_variable_mRNA.csv")
)


# 6. 初始化模型列表 ----------------------------------------------------------

models <- list()
pred_val_list <- list()
pred_test_list <- list()


# 7. Mean baseline ----------------------------------------------------------
# 因为 Y_train 已经 Z-score，所以训练集平均值约等于 0
# Mean baseline 用来证明模型不是只预测平均表达

pred_val_list$Mean_Baseline <- matrix(
  0,
  nrow = nrow(Y_val),
  ncol = ncol(Y_val),
  dimnames = list(rownames(Y_val), colnames(Y_val))
)

pred_test_list$Mean_Baseline <- matrix(
  0,
  nrow = nrow(Y_test),
  ncol = ncol(Y_test),
  dimnames = list(rownames(Y_test), colnames(Y_test))
)


# 8. Ridge ------------------------------------------------------------------

message("\n==============================")
message("Training Ridge")
message("==============================")

ridge_fit <- tryCatch(
  {
    cv.glmnet(
      x = X_train,
      y = Y_train,
      family = "mgaussian",
      alpha = 0,
      nfolds = nfolds,
      standardize = FALSE,
      type.measure = "mse"
    )
  },
  error = function(e) {
    message("Ridge 失败：", e$message)
    NULL
  }
)

if (!is.null(ridge_fit)) {

  models$Ridge <- list(
    type = "ridge",
    fit = ridge_fit
  )

  pred_val <- predict(ridge_fit, newx = X_val, s = "lambda.min")
  pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)

  pred_test <- predict(ridge_fit, newx = X_test, s = "lambda.min")
  pred_test <- to_prediction_matrix(pred_test, rownames(X_test), miRNA_names)

  pred_val_list$Ridge <- pred_val
  pred_test_list$Ridge <- pred_test
}


# 9. Shuffled baseline ------------------------------------------------------
# 打乱 mRNA-miRNA 配对关系
# 如果真实模型明显优于 shuffled baseline，说明模型确实学习了配对信号

if (run_shuffle_baseline) {

  message("\n==============================")
  message("Training Shuffled Baseline")
  message("==============================")

  set.seed(123)
  shuffle_idx <- sample(seq_len(nrow(Y_train)))
  Y_train_shuffled <- Y_train[shuffle_idx, , drop = FALSE]

  shuffled_fit <- tryCatch(
    {
      cv.glmnet(
        x = X_train,
        y = Y_train_shuffled,
        family = "mgaussian",
        alpha = 0,
        nfolds = nfolds,
        standardize = FALSE,
        type.measure = "mse"
      )
    },
    error = function(e) {
      message("Shuffled baseline 失败：", e$message)
      NULL
    }
  )

  if (!is.null(shuffled_fit)) {

    models$Shuffled_Baseline <- list(
      type = "shuffled_ridge",
      fit = shuffled_fit
    )

    pred_val <- predict(shuffled_fit, newx = X_val, s = "lambda.min")
    pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)

    pred_test <- predict(shuffled_fit, newx = X_test, s = "lambda.min")
    pred_test <- to_prediction_matrix(pred_test, rownames(X_test), miRNA_names)

    pred_val_list$Shuffled_Baseline <- pred_val
    pred_test_list$Shuffled_Baseline <- pred_test
  }
}


# 10. Elastic Net -----------------------------------------------------------

message("\n==============================")
message("Training Elastic Net")
message("==============================")

alpha_grid <- c(0.1, 0.3, 0.5, 0.7, 0.9)

enet_fits <- list()
enet_val_rmse <- c()

for (a in alpha_grid) {

  message("Elastic Net alpha = ", a)

  fit_a <- tryCatch(
    {
      cv.glmnet(
        x = X_train,
        y = Y_train,
        family = "mgaussian",
        alpha = a,
        nfolds = nfolds,
        standardize = FALSE,
        type.measure = "mse"
      )
    },
    error = function(e) {
      message("Elastic Net alpha = ", a, " 失败：", e$message)
      NULL
    }
  )

  if (!is.null(fit_a)) {

    pred_val <- predict(fit_a, newx = X_val, s = "lambda.min")
    pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)

    enet_fits[[as.character(a)]] <- fit_a
    enet_val_rmse[as.character(a)] <- sqrt(mean((Y_val - pred_val)^2, na.rm = TRUE))
  }
}

if (length(enet_fits) > 0) {

  best_alpha <- as.numeric(names(which.min(enet_val_rmse)))
  enet_fit <- enet_fits[[as.character(best_alpha)]]

  models$ElasticNet <- list(
    type = "elastic_net",
    fit = enet_fit,
    best_alpha = best_alpha,
    val_rmse_by_alpha = enet_val_rmse
  )

  pred_val <- predict(enet_fit, newx = X_val, s = "lambda.min")
  pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)

  pred_test <- predict(enet_fit, newx = X_test, s = "lambda.min")
  pred_test <- to_prediction_matrix(pred_test, rownames(X_test), miRNA_names)

  pred_val_list$ElasticNet <- pred_val
  pred_test_list$ElasticNet <- pred_test

  write_csv(
    tibble(
      alpha = as.numeric(names(enet_val_rmse)),
      validation_RMSE = as.numeric(enet_val_rmse)
    ),
    file.path(out_dir, "05_elastic_net_alpha_selection.csv")
  )
}


# 11. PLS -------------------------------------------------------------------

message("\n==============================")
message("Training PLS")
message("==============================")

pls_result <- tryCatch(
  {
    simpls_fit <- getFromNamespace("simpls.fit", "pls")

    ncomp_grid <- c(5, 10, 20, 30, 50, 80, 100)
    max_comp <- min(nrow(X_train) - 2, ncol(X_train), max(ncomp_grid))
    ncomp_grid <- ncomp_grid[ncomp_grid <= max_comp]

    if (length(ncomp_grid) == 0) {
      stop("PLS 可用成分数不足")
    }

    pls_val_rmse <- c()
    pls_coef_list <- list()

    for (k in ncomp_grid) {

      message("PLS ncomp = ", k)

      fit_k <- simpls_fit(X_train, Y_train, ncomp = k)

      coef_array <- fit_k$coefficients

      if (length(dim(coef_array)) == 3) {
        coef_k <- coef_array[, , k]
      } else {
        coef_k <- coef_array
      }

      coef_k <- as.matrix(coef_k)
      rownames(coef_k) <- colnames(X_train)
      colnames(coef_k) <- colnames(Y_train)

      pred_val <- X_val %*% coef_k
      colnames(pred_val) <- miRNA_names
      rownames(pred_val) <- rownames(X_val)

      pls_val_rmse[as.character(k)] <- sqrt(mean((Y_val - pred_val)^2, na.rm = TRUE))
      pls_coef_list[[as.character(k)]] <- coef_k
    }

    best_ncomp <- as.numeric(names(which.min(pls_val_rmse)))

    fit_final <- simpls_fit(X_train, Y_train, ncomp = best_ncomp)

    coef_array <- fit_final$coefficients

    if (length(dim(coef_array)) == 3) {
      coef_final <- coef_array[, , best_ncomp]
    } else {
      coef_final <- coef_array
    }

    coef_final <- as.matrix(coef_final)
    rownames(coef_final) <- colnames(X_train)
    colnames(coef_final) <- colnames(Y_train)

    list(
      fit = fit_final,
      coef = coef_final,
      best_ncomp = best_ncomp,
      val_rmse = pls_val_rmse
    )
  },
  error = function(e) {
    message("PLS 失败：", e$message)
    NULL
  }
)

if (!is.null(pls_result)) {

  models$PLS <- list(
    type = "pls",
    fit = pls_result$fit,
    coef = pls_result$coef,
    best_ncomp = pls_result$best_ncomp,
    val_rmse = pls_result$val_rmse
  )

  pred_val <- X_val %*% pls_result$coef
  colnames(pred_val) <- miRNA_names
  rownames(pred_val) <- rownames(X_val)

  pred_test <- X_test %*% pls_result$coef
  colnames(pred_test) <- miRNA_names
  rownames(pred_test) <- rownames(X_test)

  pred_val_list$PLS <- pred_val
  pred_test_list$PLS <- pred_test

  write_csv(
    tibble(
      ncomp = as.numeric(names(pls_result$val_rmse)),
      validation_RMSE = as.numeric(pls_result$val_rmse)
    ),
    file.path(out_dir, "06_pls_ncomp_selection.csv")
  )
}


 # 12. PCA-Ridge -------------------------------------------------------------

message("\n==============================")
message("Training PCA-Ridge")
message("==============================")

pca_ridge_result <- tryCatch(
  {
    k_grid <- c(20, 50, 100, 150, 200)
    max_k <- min(nrow(X_train) - 2, ncol(X_train), max(k_grid))
    k_grid <- k_grid[k_grid <= max_k]

    if (length(k_grid) == 0) {
      stop("PCA-Ridge 可用 PC 数不足")
    }

    pca_fit <- prcomp(X_train, center = FALSE, scale. = FALSE)

    pca_ridge_val_rmse <- c()
    pca_ridge_fits <- list()

    for (k in k_grid) {

      message("PCA-Ridge k = ", k)

      Z_train <- pca_fit$x[, 1:k, drop = FALSE]
      Z_val <- predict(pca_fit, newdata = X_val)[, 1:k, drop = FALSE]

      fit_k <- cv.glmnet(
        x = Z_train,
        y = Y_train,
        family = "mgaussian",
        alpha = 0,
        nfolds = nfolds,
        standardize = FALSE,
        type.measure = "mse"
      )

      pred_val <- predict(fit_k, newx = Z_val, s = "lambda.min")
      pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)

      pca_ridge_val_rmse[as.character(k)] <- sqrt(mean((Y_val - pred_val)^2, na.rm = TRUE))
      pca_ridge_fits[[as.character(k)]] <- fit_k
    }

    best_k <- as.numeric(names(which.min(pca_ridge_val_rmse)))

    list(
      pca = pca_fit,
      fit = pca_ridge_fits[[as.character(best_k)]],
      best_k = best_k,
      val_rmse = pca_ridge_val_rmse
    )
  },
  error = function(e) {
    message("PCA-Ridge 失败：", e$message)
    NULL
  }
)

if (!is.null(pca_ridge_result)) {

  models$PCA_Ridge <- list(
    type = "pca_ridge",
    pca = pca_ridge_result$pca,
    fit = pca_ridge_result$fit,
    best_k = pca_ridge_result$best_k,
    val_rmse = pca_ridge_result$val_rmse
  )

  Z_val <- predict(pca_ridge_result$pca, newdata = X_val)[, 1:pca_ridge_result$best_k, drop = FALSE]
  Z_test <- predict(pca_ridge_result$pca, newdata = X_test)[, 1:pca_ridge_result$best_k, drop = FALSE]

  pred_val <- predict(pca_ridge_result$fit, newx = Z_val, s = "lambda.min")
  pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)

  pred_test <- predict(pca_ridge_result$fit, newx = Z_test, s = "lambda.min")
  pred_test <- to_prediction_matrix(pred_test, rownames(X_test), miRNA_names)

  pred_val_list$PCA_Ridge <- pred_val
  pred_test_list$PCA_Ridge <- pred_test

  write_csv(
    tibble(
      k = as.numeric(names(pca_ridge_result$val_rmse)),
      validation_RMSE = as.numeric(pca_ridge_result$val_rmse)
    ),
    file.path(out_dir, "07_pca_ridge_k_selection.csv")
  )
}


# 13. PCA-ElasticNet ---------------------------------------------------------

if (run_pca_elasticnet) {

  message("\n==============================")
  message("Training PCA-ElasticNet")
  message("==============================")

  pca_elasticnet_result <- tryCatch(
    {
      k_grid <- c(20, 50, 100, 150, 200)
      alpha_grid <- c(0.1, 0.3, 0.5, 0.7, 0.9)
      max_k <- min(nrow(X_train) - 2, ncol(X_train), max(k_grid))
      k_grid <- k_grid[k_grid <= max_k]

      if (length(k_grid) == 0) {
        stop("PCA-ElasticNet 可用 PC 数不足")
      }

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

          tuning_rows[[length(tuning_rows) + 1]] <- tibble(
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

      list(
        pca = pca_fit,
        fit = best_fit,
        best_k = best_k,
        best_alpha = best_alpha,
        best_rmse = best_rmse,
        tuning = bind_rows(tuning_rows)
      )
    },
    error = function(e) {
      message("PCA-ElasticNet 失败：", e$message)
      NULL
    }
  )

  if (!is.null(pca_elasticnet_result)) {

    models$PCA_ElasticNet <- list(
      type = "pca_elasticnet",
      pca = pca_elasticnet_result$pca,
      fit = pca_elasticnet_result$fit,
      best_k = pca_elasticnet_result$best_k,
      best_alpha = pca_elasticnet_result$best_alpha,
      best_rmse = pca_elasticnet_result$best_rmse
    )

    Z_val <- predict(pca_elasticnet_result$pca, newdata = X_val)[, 1:pca_elasticnet_result$best_k, drop = FALSE]
    Z_test <- predict(pca_elasticnet_result$pca, newdata = X_test)[, 1:pca_elasticnet_result$best_k, drop = FALSE]

    pred_val <- predict(pca_elasticnet_result$fit, newx = Z_val, s = "lambda.min")
    pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)

    pred_test <- predict(pca_elasticnet_result$fit, newx = Z_test, s = "lambda.min")
    pred_test <- to_prediction_matrix(pred_test, rownames(X_test), miRNA_names)

    pred_val_list$PCA_ElasticNet <- pred_val
    pred_test_list$PCA_ElasticNet <- pred_test

    write_csv(
      pca_elasticnet_result$tuning,
      file.path(out_dir, "07b_pca_elasticnet_selection.csv")
    )
  }
}


# 14. MLP -------------------------------------------------------------------

if (run_mlp) {

  message("\n==============================")
  message("Training MLP")
  message("==============================")

  mlp_result <- tryCatch(
    {
      npc <- min(50, nrow(X_train) - 2, ncol(X_train))
      size <- 30

      if (npc < 2) {
        stop("MLP 可用 PC 数不足")
      }

      pca_mlp <- prcomp(X_train, center = FALSE, scale. = FALSE)
      Z_train <- pca_mlp$x[, 1:npc, drop = FALSE]
      Z_val <- predict(pca_mlp, newdata = X_val)[, 1:npc, drop = FALSE]
      Z_test <- predict(pca_mlp, newdata = X_test)[, 1:npc, drop = FALSE]

      mlp_fit <- nnet::nnet(
        x = Z_train,
        y = Y_train,
        size = size,
        linout = TRUE,
        skip = TRUE,
        decay = 1e-4,
        maxit = 300,
        trace = FALSE,
        MaxNWts = 500000
      )

      list(
        pca = pca_mlp,
        fit = mlp_fit,
        npc = npc,
        size = size,
        Z_val = Z_val,
        Z_test = Z_test
      )
    },
    error = function(e) {
      message("MLP 失败：", e$message)
      NULL
    }
  )

  if (!is.null(mlp_result)) {

    models$MLP <- list(
      type = "mlp",
      pca = mlp_result$pca,
      fit = mlp_result$fit,
      npc = mlp_result$npc,
      size = mlp_result$size
    )

    pred_val <- predict(mlp_result$fit, mlp_result$Z_val)
    pred_val <- as.matrix(pred_val)
    colnames(pred_val) <- miRNA_names
    rownames(pred_val) <- rownames(X_val)

    pred_test <- predict(mlp_result$fit, mlp_result$Z_test)
    pred_test <- as.matrix(pred_test)
    colnames(pred_test) <- miRNA_names
    rownames(pred_test) <- rownames(X_test)

    pred_val_list$MLP <- pred_val
    pred_test_list$MLP <- pred_test
  }
}


# 15. Pathway-Ridge ---------------------------------------------------------
# 注意：
#   需要安装 msigdbr
#   mRNA 列名最好是 gene symbol
#   如果你的 mRNA 是 Ensembl ID，这一步可能会自动跳过或报基因集匹配不足

pathway_object <- NULL

if (run_pathway_ridge) {

  message("\n==============================")
  message("Training Pathway-Ridge")
  message("==============================")

  if (!requireNamespace("msigdbr", quietly = TRUE)) {

    message("未安装 msigdbr，跳过 Pathway-Ridge。")

  } else {

    pathway_result <- tryCatch(
      {
        msig <- tryCatch(
          {
            msigdbr::msigdbr(species = "Homo sapiens", category = "H")
          },
          error = function(e1) {
            msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
          }
        )

        if (!all(c("gs_name", "gene_symbol") %in% colnames(msig))) {
          stop("msigdbr 返回格式异常")
        }

        gene_sets <- split(msig$gene_symbol, msig$gs_name)
        gene_sets <- lapply(gene_sets, unique)

        # 用训练集全部 mRNA 标准化
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

        gene_sets_used <- lapply(gene_sets, function(gs) {
          intersect(gs, colnames(X_train_all_z))
        })

        gene_sets_used <- gene_sets_used[
          lengths(gene_sets_used) >= 5 &
            lengths(gene_sets_used) <= 500
        ]

        if (length(gene_sets_used) == 0) {
          stop("没有可用 pathway，请检查 mRNA 列名是否为 gene symbol")
        }

        P_train_raw <- matrix(
          NA_real_,
          nrow = nrow(X_train_all_z),
          ncol = length(gene_sets_used),
          dimnames = list(rownames(X_train_all_z), names(gene_sets_used))
        )

        P_val_raw <- matrix(
          NA_real_,
          nrow = nrow(X_val_all_z),
          ncol = length(gene_sets_used),
          dimnames = list(rownames(X_val_all_z), names(gene_sets_used))
        )

        P_test_raw <- matrix(
          NA_real_,
          nrow = nrow(X_test_all_z),
          ncol = length(gene_sets_used),
          dimnames = list(rownames(X_test_all_z), names(gene_sets_used))
        )

        for (i in seq_along(gene_sets_used)) {

          gs <- gene_sets_used[[i]]

          P_train_raw[, i] <- rowMeans(X_train_all_z[, gs, drop = FALSE], na.rm = TRUE)
          P_val_raw[, i] <- rowMeans(X_val_all_z[, gs, drop = FALSE], na.rm = TRUE)
          P_test_raw[, i] <- rowMeans(X_test_all_z[, gs, drop = FALSE], na.rm = TRUE)
        }

        P_train_raw[!is.finite(P_train_raw)] <- 0
        P_val_raw[!is.finite(P_val_raw)] <- 0
        P_test_raw[!is.finite(P_test_raw)] <- 0

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

        pathway_fit <- cv.glmnet(
          x = P_train,
          y = Y_train,
          family = "mgaussian",
          alpha = 0,
          nfolds = nfolds,
          standardize = FALSE,
          type.measure = "mse"
        )

        list(
          fit = pathway_fit,
          gene_sets = gene_sets,
          gene_sets_used = gene_sets_used,
          x_all_center = x_all_center,
          x_all_scale = x_all_scale,
          pathway_center = pathway_center,
          pathway_scale = pathway_scale,
          pathway_names = colnames(P_train),
          P_val = P_val,
          P_test = P_test
        )
      },
      error = function(e) {
        message("Pathway-Ridge 失败：", e$message)
        NULL
      }
    )

    if (!is.null(pathway_result)) {

      models$Pathway_Ridge <- list(
        type = "pathway_ridge",
        fit = pathway_result$fit
      )

      pathway_object <- list(
        gene_sets = pathway_result$gene_sets,
        gene_sets_used = pathway_result$gene_sets_used,
        x_all_center = pathway_result$x_all_center,
        x_all_scale = pathway_result$x_all_scale,
        pathway_center = pathway_result$pathway_center,
        pathway_scale = pathway_result$pathway_scale,
        pathway_names = pathway_result$pathway_names
      )

      pred_val <- predict(pathway_result$fit, newx = pathway_result$P_val, s = "lambda.min")
      pred_val <- to_prediction_matrix(pred_val, rownames(X_val), miRNA_names)

      pred_test <- predict(pathway_result$fit, newx = pathway_result$P_test, s = "lambda.min")
      pred_test <- to_prediction_matrix(pred_test, rownames(X_test), miRNA_names)

      pred_val_list$Pathway_Ridge <- pred_val
      pred_test_list$Pathway_Ridge <- pred_test

      # Pathway coefficient interpretation
      coef_list <- coef(pathway_result$fit, s = "lambda.min")

      pathway_coef <- sapply(coef_list, function(x) {
        as.numeric(x[-1, 1])
      })

      pathway_coef <- as.matrix(pathway_coef)
      rownames(pathway_coef) <- pathway_result$pathway_names
      colnames(pathway_coef) <- miRNA_names

      pathway_summary <- tibble(
        pathway = rownames(pathway_coef),
        mean_abs_coef = rowMeans(abs(pathway_coef), na.rm = TRUE),
        max_abs_coef = apply(abs(pathway_coef), 1, max, na.rm = TRUE)
      ) %>%
        arrange(desc(mean_abs_coef))

      write_csv(
        pathway_summary,
        file.path(out_dir, "08_pathway_ridge_interpretation_summary.csv")
      )

      write_csv(
        data.frame(
          pathway = rownames(pathway_coef),
          pathway_coef,
          check.names = FALSE
        ),
        file.path(out_dir, "09_pathway_ridge_coef_matrix.csv")
      )
    }
  }
}


# 15. validation 和 test 模型评估 ------------------------------------------

if (length(pred_val_list) == 0 || length(pred_test_list) == 0) {
  stop("没有任何模型成功，请检查数据格式。")
}

metric_val <- list()

for (model_name in names(pred_val_list)) {

  real <- Y_val
  pred <- pred_val_list[[model_name]]

  rss <- colSums((real - pred)^2, na.rm = TRUE)
  tss <- colSums(
    sweep(real, 2, colMeans(real, na.rm = TRUE), "-")^2,
    na.rm = TRUE
  )

  R2 <- 1 - rss / tss
  R2[!is.finite(R2)] <- NA_real_

  Pearson <- rep(NA_real_, ncol(real))

  for (j in seq_len(ncol(real))) {
    if (sd(real[, j], na.rm = TRUE) > 1e-8 &&
        sd(pred[, j], na.rm = TRUE) > 1e-8) {
      Pearson[j] <- suppressWarnings(
        cor(real[, j], pred[, j], use = "complete.obs")
      )
    }
  }

  RMSE <- sqrt(colMeans((real - pred)^2, na.rm = TRUE))

  metric_val[[model_name]] <- tibble(
    dataset = "validation",
    model = model_name,
    miRNA = colnames(real),
    R2 = R2,
    Pearson = Pearson,
    RMSE = RMSE
  )
}

metric_val <- bind_rows(metric_val)

metric_test <- list()

for (model_name in names(pred_test_list)) {

  real <- Y_test
  pred <- pred_test_list[[model_name]]

  rss <- colSums((real - pred)^2, na.rm = TRUE)
  tss <- colSums(
    sweep(real, 2, colMeans(real, na.rm = TRUE), "-")^2,
    na.rm = TRUE
  )

  R2 <- 1 - rss / tss
  R2[!is.finite(R2)] <- NA_real_

  Pearson <- rep(NA_real_, ncol(real))

  for (j in seq_len(ncol(real))) {
    if (sd(real[, j], na.rm = TRUE) > 1e-8 &&
        sd(pred[, j], na.rm = TRUE) > 1e-8) {
      Pearson[j] <- suppressWarnings(
        cor(real[, j], pred[, j], use = "complete.obs")
      )
    }
  }

  RMSE <- sqrt(colMeans((real - pred)^2, na.rm = TRUE))

  metric_test[[model_name]] <- tibble(
    dataset = "test",
    model = model_name,
    miRNA = colnames(real),
    R2 = R2,
    Pearson = Pearson,
    RMSE = RMSE
  )
}

metric_test <- bind_rows(metric_test)

global_val <- metric_val %>%
  group_by(model) %>%
  summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    median_R2 = median(R2, na.rm = TRUE),
    mean_Pearson = mean(Pearson, na.rm = TRUE),
    median_Pearson = median(Pearson, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    median_RMSE = median(RMSE, na.rm = TRUE),
    n_miRNA = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_Pearson), desc(mean_R2), mean_RMSE)

global_test <- metric_test %>%
  group_by(model) %>%
  summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    median_R2 = median(R2, na.rm = TRUE),
    mean_Pearson = mean(Pearson, na.rm = TRUE),
    median_Pearson = median(Pearson, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    median_RMSE = median(RMSE, na.rm = TRUE),
    n_miRNA = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_Pearson), desc(mean_R2), mean_RMSE)

# 只用 validation 选择最佳模型
valid_candidate <- intersect(candidate_models, global_val$model)

if (length(valid_candidate) == 0) {
  stop("候选模型均未成功训练。")
}

best_model <- global_val %>%
  filter(model %in% valid_candidate) %>%
  arrange(desc(median_Pearson), desc(mean_R2), mean_RMSE) %>%
  slice(1) %>%
  pull(model)

message("\n根据 validation 选择的最佳模型：", best_model)

write_csv(metric_val, file.path(out_dir, "10_validation_metrics_by_miRNA.csv"))
write_csv(global_val, file.path(out_dir, "11_validation_global_metrics.csv"))

write_csv(metric_test, file.path(out_dir, "12_test_metrics_by_miRNA.csv"))
write_csv(global_test, file.path(out_dir, "13_test_global_metrics.csv"))


# 16. 输出最佳模型 virtual miRNA --------------------------------------------

best_pred_val_z <- pred_val_list[[best_model]]
best_pred_test_z <- pred_test_list[[best_model]]

best_pred_val_expr <- sweep(best_pred_val_z, 2, y_scale, "*")
best_pred_val_expr <- sweep(best_pred_val_expr, 2, y_center, "+")

best_pred_test_expr <- sweep(best_pred_test_z, 2, y_scale, "*")
best_pred_test_expr <- sweep(best_pred_test_expr, 2, y_center, "+")

write_csv(
  data.frame(
    patient_id = rownames(best_pred_val_z),
    best_pred_val_z,
    check.names = FALSE
  ),
  file.path(out_dir, paste0("14_validation_virtual_miRNA_z_", best_model, ".csv"))
)

write_csv(
  data.frame(
    patient_id = rownames(best_pred_test_z),
    best_pred_test_z,
    check.names = FALSE
  ),
  file.path(out_dir, paste0("15_test_virtual_miRNA_z_", best_model, ".csv"))
)

write_csv(
  data.frame(
    patient_id = rownames(best_pred_test_expr),
    best_pred_test_expr,
    check.names = FALSE
  ),
  file.path(out_dir, paste0("16_test_virtual_miRNA_expression_", best_model, ".csv"))
)


# 17. miRNA-level reliability -----------------------------------------------
# 思路：
#   在 validation 上估计每个 miRNA 的可靠性
#   然后在 test 上检验 reliability 是否与 error 负相关

miRNA_reliability <- metric_val %>%
  filter(model == best_model) %>%
  mutate(
    R2_clamped = pmin(pmax(R2, 0), 1),
    Pearson_clamped = pmin(pmax(Pearson, 0), 1),
    miRNA_level_reliability = 0.5 * R2_clamped + 0.5 * Pearson_clamped,
    validation_prediction_error = RMSE
  ) %>%
  select(
    miRNA,
    R2,
    Pearson,
    RMSE,
    miRNA_level_reliability,
    validation_prediction_error
  ) %>%
  arrange(desc(miRNA_level_reliability))

miRNA_test_error <- metric_test %>%
  filter(model == best_model) %>%
  select(
    miRNA,
    test_R2 = R2,
    test_Pearson = Pearson,
    test_RMSE = RMSE
  )

miRNA_reliability_test_check <- miRNA_reliability %>%
  left_join(miRNA_test_error, by = "miRNA")

miRNA_rel_test_error_cor <- suppressWarnings(
  cor(
    miRNA_reliability_test_check$miRNA_level_reliability,
    miRNA_reliability_test_check$test_RMSE,
    method = "spearman",
    use = "complete.obs"
  )
)

write_csv(
  miRNA_reliability,
  file.path(out_dir, "17_miRNA_level_reliability_from_validation.csv")
)

write_csv(
  miRNA_reliability_test_check,
  file.path(out_dir, "18_miRNA_reliability_test_error_check.csv")
)


# 18. sample-level reliability ----------------------------------------------
# 两部分：
#   distance_score：样本是否接近训练集分布
#   agreement_score：多个模型预测是否一致
#
# 注意：
#   agreement 的 min / max 来自 validation
#   test 和未来 private_mRNA 都复用 validation 参数

n_pc_distance <- min(20, nrow(X_train) - 2, ncol(X_train))

pca_distance <- prcomp(X_train, center = FALSE, scale. = FALSE)

train_pc <- pca_distance$x[, 1:n_pc_distance, drop = FALSE]
val_pc <- predict(pca_distance, newdata = X_val)[, 1:n_pc_distance, drop = FALSE]
test_pc <- predict(pca_distance, newdata = X_test)[, 1:n_pc_distance, drop = FALSE]

distance_center <- colMeans(train_pc)
distance_variance <- apply(train_pc, 2, var)
distance_variance[!is.finite(distance_variance) | distance_variance < 1e-8] <- 1

d_train <- rowSums(
  sweep(train_pc, 2, distance_center, "-")^2 /
    matrix(distance_variance, nrow(train_pc), n_pc_distance, byrow = TRUE)
)

d_val <- rowSums(
  sweep(val_pc, 2, distance_center, "-")^2 /
    matrix(distance_variance, nrow(val_pc), n_pc_distance, byrow = TRUE)
)

d_test <- rowSums(
  sweep(test_pc, 2, distance_center, "-")^2 /
    matrix(distance_variance, nrow(test_pc), n_pc_distance, byrow = TRUE)
)

distance_q95 <- as.numeric(quantile(d_train, 0.95, na.rm = TRUE))

if (!is.finite(distance_q95) || distance_q95 <= 0) {
  distance_q95 <- median(d_train, na.rm = TRUE)
}

distance_score_val <- exp(-d_val / distance_q95)
distance_score_test <- exp(-d_test / distance_q95)

available_prediction_models <- intersect(candidate_models, names(pred_val_list))

if (length(available_prediction_models) >= 2) {

  pred_array_val <- array(
    NA_real_,
    dim = c(nrow(Y_val), ncol(Y_val), length(available_prediction_models)),
    dimnames = list(
      rownames(Y_val),
      colnames(Y_val),
      available_prediction_models
    )
  )

  pred_array_test <- array(
    NA_real_,
    dim = c(nrow(Y_test), ncol(Y_test), length(available_prediction_models)),
    dimnames = list(
      rownames(Y_test),
      colnames(Y_test),
      available_prediction_models
    )
  )

  for (i in seq_along(available_prediction_models)) {
    m <- available_prediction_models[i]
    pred_array_val[, , i] <- pred_val_list[[m]]
    pred_array_test[, , i] <- pred_test_list[[m]]
  }

  model_disagreement_val <- apply(pred_array_val, 1, function(x) {
    mean(apply(x, 1, sd, na.rm = TRUE), na.rm = TRUE)
  })

  model_disagreement_test <- apply(pred_array_test, 1, function(x) {
    mean(apply(x, 1, sd, na.rm = TRUE), na.rm = TRUE)
  })

  disagreement_min <- min(model_disagreement_val, na.rm = TRUE)
  disagreement_max <- max(model_disagreement_val, na.rm = TRUE)

  if (!is.finite(disagreement_min) ||
      !is.finite(disagreement_max) ||
      abs(disagreement_max - disagreement_min) < 1e-12) {

    agreement_score_val <- rep(NA_real_, length(model_disagreement_val))
    agreement_score_test <- rep(NA_real_, length(model_disagreement_test))

  } else {

    agreement_score_val <- 1 - (model_disagreement_val - disagreement_min) /
      (disagreement_max - disagreement_min)

    agreement_score_test <- 1 - (model_disagreement_test - disagreement_min) /
      (disagreement_max - disagreement_min)

    agreement_score_val <- pmin(pmax(agreement_score_val, 0), 1)
    agreement_score_test <- pmin(pmax(agreement_score_test, 0), 1)
  }

} else {

  model_disagreement_val <- rep(NA_real_, nrow(Y_val))
  model_disagreement_test <- rep(NA_real_, nrow(Y_test))
  agreement_score_val <- rep(NA_real_, nrow(Y_val))
  agreement_score_test <- rep(NA_real_, nrow(Y_test))

  disagreement_min <- NA_real_
  disagreement_max <- NA_real_
}

sample_error_val <- sqrt(rowMeans((Y_val - best_pred_val_z)^2, na.rm = TRUE))
sample_error_test <- sqrt(rowMeans((Y_test - best_pred_test_z)^2, na.rm = TRUE))

sample_reliability_val <- tibble(
  patient_id = rownames(Y_val),
  dataset = "validation",
  distance_score = as.numeric(distance_score_val),
  model_disagreement = as.numeric(model_disagreement_val),
  agreement_score = as.numeric(agreement_score_val),
  sample_level_reliability = rowMeans(
    cbind(distance_score_val, agreement_score_val),
    na.rm = TRUE
  ),
  sample_prediction_error = sample_error_val
)

sample_reliability_test <- tibble(
  patient_id = rownames(Y_test),
  dataset = "test",
  distance_score = as.numeric(distance_score_test),
  model_disagreement = as.numeric(model_disagreement_test),
  agreement_score = as.numeric(agreement_score_test),
  sample_level_reliability = rowMeans(
    cbind(distance_score_test, agreement_score_test),
    na.rm = TRUE
  ),
  sample_prediction_error = sample_error_test
)

sample_rel_val_error_cor <- suppressWarnings(
  cor(
    sample_reliability_val$sample_level_reliability,
    sample_reliability_val$sample_prediction_error,
    method = "spearman",
    use = "complete.obs"
  )
)

sample_rel_test_error_cor <- suppressWarnings(
  cor(
    sample_reliability_test$sample_level_reliability,
    sample_reliability_test$sample_prediction_error,
    method = "spearman",
    use = "complete.obs"
  )
)

write_csv(
  sample_reliability_val,
  file.path(out_dir, "19_sample_level_reliability_validation.csv")
)

write_csv(
  sample_reliability_test,
  file.path(out_dir, "20_sample_level_reliability_test.csv")
)

write_csv(
  tibble(
    item = c(
      "miRNA_reliability_vs_test_error_spearman",
      "sample_reliability_vs_validation_error_spearman",
      "sample_reliability_vs_test_error_spearman",
      "agreement_disagreement_min_from_validation",
      "agreement_disagreement_max_from_validation"
    ),
    value = c(
      miRNA_rel_test_error_cor,
      sample_rel_val_error_cor,
      sample_rel_test_error_cor,
      disagreement_min,
      disagreement_max
    )
  ),
  file.path(out_dir, "21_reliability_calibration_summary.csv")
)



# 19. Downstream preservation analysis by miRNA reliability -----------------
# 说明：
#   原来的 downstream preservation 是基于全部 miRNA。
#   这里改成四套结果：
#     All
#     Rel >= 0.5
#     Rel >= 0.6
#     Rel >= 0.7
#   这样可以证明：miRNA-level reliability 越高，virtual miRNA 越能保留下游信号。

real_z <- Y_test
virtual_z <- best_pred_test_z

real_expr <- Y_test_raw
virtual_expr <- best_pred_test_expr

rel_out_dir <- file.path(out_dir, "downstream_by_miRNA_reliability")
rel_fig_dir <- file.path(fig_dir, "downstream_by_miRNA_reliability")

dir.create(rel_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rel_fig_dir, showWarnings = FALSE, recursive = TRUE)

make_miRNA_set <- function(miRNA_reliability, all_miRNA, cutoff = NA_real_) {
  if (is.na(cutoff)) {
    return(all_miRNA)
  }

  miRNA_reliability %>%
    filter(miRNA_level_reliability >= cutoff) %>%
    pull(miRNA) %>%
    intersect(all_miRNA)
}

miRNA_sets <- list(
  All = make_miRNA_set(miRNA_reliability, colnames(real_z), NA_real_),
  Rel_ge_0_5 = make_miRNA_set(miRNA_reliability, colnames(real_z), 0.5),
  Rel_ge_0_6 = make_miRNA_set(miRNA_reliability, colnames(real_z), 0.6),
  Rel_ge_0_7 = make_miRNA_set(miRNA_reliability, colnames(real_z), 0.7)
)

miRNA_set_summary <- tibble(
  miRNA_set = names(miRNA_sets),
  cutoff = c(NA_real_, 0.5, 0.6, 0.7),
  n_miRNA = as.integer(sapply(miRNA_sets, length))
)

write_csv(
  miRNA_set_summary,
  file.path(rel_out_dir, "00_miRNA_reliability_threshold_summary.csv")
)

# ARI / NMI 小函数 ----------------------------------------------------------

ari_score <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)

  comb2 <- function(z) z * (z - 1) / 2

  if (n < 2) return(NA_real_)

  sum_nij <- sum(comb2(tab))
  sum_ai <- sum(comb2(rowSums(tab)))
  sum_bj <- sum(comb2(colSums(tab)))

  expected <- sum_ai * sum_bj / comb2(n)
  max_index <- (sum_ai + sum_bj) / 2

  if (abs(max_index - expected) < 1e-12) return(NA_real_)

  (sum_nij - expected) / (max_index - expected)
}

nmi_score <- function(x, y) {
  tab <- table(x, y)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)

  nz <- pxy > 0

  if (!any(nz)) return(NA_real_)

  r <- row(pxy)[nz]
  c <- col(pxy)[nz]

  mi <- sum(pxy[nz] * log(pxy[nz] / (px[r] * py[c])))

  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))

  if (hx < 1e-12 || hy < 1e-12) return(NA_real_)

  mi / sqrt(hx * hy)
}

get_downstream_group <- function(clinical_test) {
  diff_group <- NULL
  group_source <- NA_character_

  if ("stage" %in% colnames(clinical_test)) {
    st <- toupper(as.character(clinical_test$stage))
    st <- gsub("STAGE", "", st)
    st <- gsub("[[:space:]]", "", st)

    group_tmp <- rep(NA_character_, length(st))

    early_stage <- c("I", "IA", "IB", "IC", "II", "IIA", "IIB", "IIC")
    late_stage <- c("III", "IIIA", "IIIB", "IIIC", "IV", "IVA", "IVB")

    group_tmp[st %in% early_stage] <- "Early"
    group_tmp[st %in% late_stage] <- "Late"

    if (sum(group_tmp == "Early", na.rm = TRUE) >= 5 &&
        sum(group_tmp == "Late", na.rm = TRUE) >= 5) {
      diff_group <- group_tmp
      group_source <- "stage_Early_vs_Late"
    }
  }

  if (is.null(diff_group) && "OS_status" %in% colnames(clinical_test)) {
    status_tmp <- as.character(clinical_test$OS_status)

    group_tmp <- ifelse(
      status_tmp %in% c("1", "Dead", "DEAD", "death", "Death", "event", "Event"),
      "Event",
      "No_event"
    )

    if (sum(group_tmp == "Event", na.rm = TRUE) >= 5 &&
        sum(group_tmp == "No_event", na.rm = TRUE) >= 5) {
      diff_group <- group_tmp
      group_source <- "OS_status_Event_vs_No_event"
    }
  }

  list(group = diff_group, source = group_source)
}

run_downstream_for_miRNA_set <- function(
    set_name,
    keep_miRNA,
    real_z,
    virtual_z,
    real_expr,
    virtual_expr,
    clinical_test,
    best_model,
    rel_out_dir,
    rel_fig_dir
) {

  message("\n==============================")
  message("Downstream analysis: ", set_name)
  message("miRNA number: ", length(keep_miRNA))
  message("==============================")

  if (length(keep_miRNA) < 2) {
    message("miRNA 数量少于 2，跳过该阈值：", set_name)
    return(list(
      summary = tibble(
        miRNA_set = set_name,
        n_miRNA = length(keep_miRNA),
        median_sample_correlation = NA_real_,
        mean_sample_correlation = NA_real_,
        correlation_matrix_preservation = NA_real_,
        differential_logFC_correlation = NA_real_,
        differential_top50_overlap = NA_real_,
        cox_logHR_correlation = NA_real_
      ),
      subtype = NULL
    ))
  }

  prefix <- paste0(set_name, "_", best_model)

  real_z_sub <- real_z[, keep_miRNA, drop = FALSE]
  virtual_z_sub <- virtual_z[, keep_miRNA, drop = FALSE]

  real_expr_sub <- real_expr[, keep_miRNA, drop = FALSE]
  virtual_expr_sub <- virtual_expr[, keep_miRNA, drop = FALSE]

  # 19.1 PCA preservation --------------------------------------------------
  combined_mat <- rbind(real_z_sub, virtual_z_sub)
  pca_combined <- prcomp(combined_mat, center = TRUE, scale. = FALSE)

  pca_df <- tibble(
    patient_id = rep(rownames(real_z_sub), 2),
    source = rep(c("Real", "Virtual"), each = nrow(real_z_sub)),
    PC1 = pca_combined$x[, 1],
    PC2 = pca_combined$x[, 2],
    miRNA_set = set_name,
    n_miRNA = length(keep_miRNA)
  )

  write_csv(
    pca_df,
    file.path(rel_out_dir, paste0(prefix, "_01_pca_real_vs_virtual.csv"))
  )

  p_pca <- ggplot(pca_df, aes(PC1, PC2, shape = source)) +
    geom_point(alpha = 0.8) +
    theme_bw() +
    labs(
      title = paste0("PCA preservation: ", set_name),
      subtitle = paste0("Best model = ", best_model, "; n miRNA = ", length(keep_miRNA)),
      x = "PC1",
      y = "PC2"
    )

  ggsave(
    file.path(rel_fig_dir, paste0(prefix, "_01_pca_real_vs_virtual.pdf")),
    p_pca,
    width = 5.5,
    height = 4.5
  )

  # 19.2 sample-wise real-virtual correlation -----------------------------
  sample_cor <- rep(NA_real_, nrow(real_z_sub))

  for (i in seq_len(nrow(real_z_sub))) {
    if (sd(real_z_sub[i, ], na.rm = TRUE) > 1e-8 &&
        sd(virtual_z_sub[i, ], na.rm = TRUE) > 1e-8) {
      sample_cor[i] <- suppressWarnings(
        cor(real_z_sub[i, ], virtual_z_sub[i, ], use = "complete.obs")
      )
    }
  }

  sample_cor_df <- tibble(
    patient_id = rownames(real_z_sub),
    miRNA_set = set_name,
    n_miRNA = length(keep_miRNA),
    real_virtual_sample_correlation = sample_cor
  )

  write_csv(
    sample_cor_df,
    file.path(rel_out_dir, paste0(prefix, "_02_sample_correlation.csv"))
  )

  p_sample_cor <- ggplot(sample_cor_df, aes(x = real_virtual_sample_correlation)) +
    geom_histogram(bins = 30, alpha = 0.8) +
    theme_bw() +
    labs(
      title = paste0("Sample-wise real-virtual correlation: ", set_name),
      subtitle = paste0("Median r = ", round(median(sample_cor, na.rm = TRUE), 3)),
      x = "Sample-wise correlation",
      y = "Number of samples"
    )

  ggsave(
    file.path(rel_fig_dir, paste0(prefix, "_02_sample_correlation_histogram.pdf")),
    p_sample_cor,
    width = 5,
    height = 4
  )

  # 19.3 miRNA-miRNA correlation preservation -----------------------------
  cor_matrix_preservation <- NA_real_

  if (ncol(real_z_sub) >= 3) {
    cor_real <- cor(real_z_sub, use = "pairwise.complete.obs")
    cor_virtual <- cor(virtual_z_sub, use = "pairwise.complete.obs")

    idx_cor <- upper.tri(cor_real)

    cor_matrix_preservation <- suppressWarnings(
      cor(
        cor_real[idx_cor],
        cor_virtual[idx_cor],
        use = "complete.obs"
      )
    )

    cor_plot_df <- tibble(
      real_correlation = cor_real[idx_cor],
      virtual_correlation = cor_virtual[idx_cor]
    ) %>%
      filter(is.finite(real_correlation), is.finite(virtual_correlation))

    if (nrow(cor_plot_df) > 50000) {
      set.seed(123)
      cor_plot_df <- cor_plot_df %>% slice_sample(n = 50000)
    }

    write_csv(
      tibble(
        miRNA_set = set_name,
        n_miRNA = length(keep_miRNA),
        item = "correlation_matrix_preservation",
        value = cor_matrix_preservation
      ),
      file.path(rel_out_dir, paste0(prefix, "_03_correlation_matrix_preservation_summary.csv"))
    )

    write_csv(
      cor_plot_df,
      file.path(rel_out_dir, paste0(prefix, "_03_correlation_matrix_preservation_plot_data.csv"))
    )

    p_cor <- ggplot(cor_plot_df, aes(real_correlation, virtual_correlation)) +
      geom_point(alpha = 0.2) +
      geom_smooth(method = "lm", se = FALSE) +
      theme_bw() +
      labs(
        title = paste0("miRNA-miRNA correlation preservation: ", set_name),
        subtitle = paste0("r = ", round(cor_matrix_preservation, 3),
                          "; n miRNA = ", length(keep_miRNA)),
        x = "Real miRNA correlation",
        y = "Virtual miRNA correlation"
      )

    ggsave(
      file.path(rel_fig_dir, paste0(prefix, "_03_correlation_matrix_preservation.pdf")),
      p_cor,
      width = 5,
      height = 4
    )
  }

  # 19.4 differential logFC preservation ----------------------------------
  diff_logfc_cor <- NA_real_
  top_overlap <- NA_real_

  group_info <- get_downstream_group(clinical_test)
  diff_group <- group_info$group

  if (!is.null(diff_group)) {
    group_factor <- factor(diff_group)
    keep_group <- !is.na(group_factor)
    group_factor <- droplevels(group_factor[keep_group])

    if (nlevels(group_factor) == 2) {
      real_mat_diff <- real_expr_sub[keep_group, , drop = FALSE]
      virtual_mat_diff <- virtual_expr_sub[keep_group, , drop = FALSE]

      g1 <- levels(group_factor)[1]
      g2 <- levels(group_factor)[2]

      real_logFC <- colMeans(real_mat_diff[group_factor == g2, , drop = FALSE], na.rm = TRUE) -
        colMeans(real_mat_diff[group_factor == g1, , drop = FALSE], na.rm = TRUE)

      virtual_logFC <- colMeans(virtual_mat_diff[group_factor == g2, , drop = FALSE], na.rm = TRUE) -
        colMeans(virtual_mat_diff[group_factor == g1, , drop = FALSE], na.rm = TRUE)

      real_pvalue <- rep(NA_real_, ncol(real_mat_diff))
      virtual_pvalue <- rep(NA_real_, ncol(virtual_mat_diff))

      for (j in seq_len(ncol(real_mat_diff))) {
        real_pvalue[j] <- tryCatch(
          {
            t.test(
              real_mat_diff[group_factor == g2, j],
              real_mat_diff[group_factor == g1, j]
            )$p.value
          },
          error = function(e) NA_real_
        )

        virtual_pvalue[j] <- tryCatch(
          {
            t.test(
              virtual_mat_diff[group_factor == g2, j],
              virtual_mat_diff[group_factor == g1, j]
            )$p.value
          },
          error = function(e) NA_real_
        )
      }

      diff_compare <- tibble(
        miRNA_set = set_name,
        n_miRNA = length(keep_miRNA),
        group_source = group_info$source,
        miRNA = colnames(real_mat_diff),
        group1 = g1,
        group2 = g2,
        real_logFC = real_logFC,
        virtual_logFC = virtual_logFC,
        real_pvalue = real_pvalue,
        virtual_pvalue = virtual_pvalue,
        real_FDR = p.adjust(real_pvalue, method = "BH"),
        virtual_FDR = p.adjust(virtual_pvalue, method = "BH"),
        abs_logFC_error = abs(real_logFC - virtual_logFC)
      )

      diff_logfc_cor <- suppressWarnings(
        cor(
          diff_compare$real_logFC,
          diff_compare$virtual_logFC,
          use = "complete.obs"
        )
      )

      top_n <- min(50, nrow(diff_compare))

      if (top_n > 0) {
        top_real <- diff_compare %>%
          arrange(real_pvalue) %>%
          slice_head(n = top_n) %>%
          pull(miRNA)

        top_virtual <- diff_compare %>%
          arrange(virtual_pvalue) %>%
          slice_head(n = top_n) %>%
          pull(miRNA)

        top_overlap <- length(intersect(top_real, top_virtual)) / top_n
      }

      write_csv(
        diff_compare,
        file.path(rel_out_dir, paste0(prefix, "_04_differential_logFC_preservation.csv"))
      )

      p_diff <- ggplot(diff_compare, aes(real_logFC, virtual_logFC)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE) +
        theme_bw() +
        labs(
          title = paste0("Differential logFC preservation: ", set_name),
          subtitle = paste0("r = ", round(diff_logfc_cor, 3),
                            "; top50 overlap = ", round(top_overlap, 3)),
          x = "Real miRNA logFC",
          y = "Virtual miRNA logFC"
        )

      ggsave(
        file.path(rel_fig_dir, paste0(prefix, "_04_differential_logFC_preservation.pdf")),
        p_diff,
        width = 5,
        height = 4
      )
    }
  }

  # 19.5 Cox HR preservation ----------------------------------------------
  cox_cor <- NA_real_

  if (all(c("OS_time", "OS_status") %in% colnames(clinical_test))) {
    time <- suppressWarnings(as.numeric(clinical_test$OS_time))
    status_raw <- as.character(clinical_test$OS_status)

    status <- rep(NA_real_, length(status_raw))
    status[status_raw %in% c("1", "Dead", "DEAD", "death", "Death", "event", "Event")] <- 1
    status[status_raw %in% c("0", "Alive", "ALIVE", "living", "Living", "No_event", "no_event")] <- 0

    status_num <- suppressWarnings(as.numeric(status_raw))
    status[is.na(status) & status_num %in% c(0, 1)] <- status_num[is.na(status) & status_num %in% c(0, 1)]

    real_center_cox <- colMeans(real_expr_sub, na.rm = TRUE)
    real_scale_cox <- apply(real_expr_sub, 2, sd, na.rm = TRUE)
    real_scale_cox[!is.finite(real_scale_cox) | real_scale_cox < 1e-8] <- 1

    virtual_center_cox <- colMeans(virtual_expr_sub, na.rm = TRUE)
    virtual_scale_cox <- apply(virtual_expr_sub, 2, sd, na.rm = TRUE)
    virtual_scale_cox[!is.finite(virtual_scale_cox) | virtual_scale_cox < 1e-8] <- 1

    real_mat_cox_z <- sweep(real_expr_sub, 2, real_center_cox, "-")
    real_mat_cox_z <- sweep(real_mat_cox_z, 2, real_scale_cox, "/")

    virtual_mat_cox_z <- sweep(virtual_expr_sub, 2, virtual_center_cox, "-")
    virtual_mat_cox_z <- sweep(virtual_mat_cox_z, 2, virtual_scale_cox, "/")

    cox_real_res <- list()
    cox_virtual_res <- list()

    for (j in seq_len(ncol(real_mat_cox_z))) {
      mir <- colnames(real_mat_cox_z)[j]

      df_real <- tibble(
        time = time,
        status = status,
        expr = real_mat_cox_z[, j]
      ) %>%
        filter(
          is.finite(time),
          is.finite(status),
          is.finite(expr),
          time > 0
        )

      if (nrow(df_real) >= 20 && length(unique(df_real$status)) >= 2) {
        fit_real <- tryCatch(
          coxph(Surv(time, status) ~ expr, data = df_real),
          error = function(e) NULL
        )

        if (!is.null(fit_real)) {
          s <- summary(fit_real)
          cox_real_res[[mir]] <- tibble(
            miRNA = mir,
            real_HR = s$coefficients[1, "exp(coef)"],
            real_logHR = s$coefficients[1, "coef"],
            real_pvalue = s$coefficients[1, "Pr(>|z|)"]
          )
        } else {
          cox_real_res[[mir]] <- tibble(
            miRNA = mir,
            real_HR = NA_real_,
            real_logHR = NA_real_,
            real_pvalue = NA_real_
          )
        }
      } else {
        cox_real_res[[mir]] <- tibble(
          miRNA = mir,
          real_HR = NA_real_,
          real_logHR = NA_real_,
          real_pvalue = NA_real_
        )
      }

      df_virtual <- tibble(
        time = time,
        status = status,
        expr = virtual_mat_cox_z[, j]
      ) %>%
        filter(
          is.finite(time),
          is.finite(status),
          is.finite(expr),
          time > 0
        )

      if (nrow(df_virtual) >= 20 && length(unique(df_virtual$status)) >= 2) {
        fit_virtual <- tryCatch(
          coxph(Surv(time, status) ~ expr, data = df_virtual),
          error = function(e) NULL
        )

        if (!is.null(fit_virtual)) {
          s <- summary(fit_virtual)
          cox_virtual_res[[mir]] <- tibble(
            miRNA = mir,
            virtual_HR = s$coefficients[1, "exp(coef)"],
            virtual_logHR = s$coefficients[1, "coef"],
            virtual_pvalue = s$coefficients[1, "Pr(>|z|)"]
          )
        } else {
          cox_virtual_res[[mir]] <- tibble(
            miRNA = mir,
            virtual_HR = NA_real_,
            virtual_logHR = NA_real_,
            virtual_pvalue = NA_real_
          )
        }
      } else {
        cox_virtual_res[[mir]] <- tibble(
          miRNA = mir,
          virtual_HR = NA_real_,
          virtual_logHR = NA_real_,
          virtual_pvalue = NA_real_
        )
      }
    }

    cox_real_df <- bind_rows(cox_real_res) %>%
      mutate(real_FDR = p.adjust(real_pvalue, method = "BH"))

    cox_virtual_df <- bind_rows(cox_virtual_res) %>%
      mutate(virtual_FDR = p.adjust(virtual_pvalue, method = "BH"))

    cox_compare <- cox_real_df %>%
      left_join(cox_virtual_df, by = "miRNA") %>%
      mutate(
        miRNA_set = set_name,
        n_miRNA = length(keep_miRNA)
      ) %>%
      select(miRNA_set, n_miRNA, everything())

    cox_cor <- suppressWarnings(
      cor(
        cox_compare$real_logHR,
        cox_compare$virtual_logHR,
        use = "complete.obs"
      )
    )

    write_csv(
      cox_compare,
      file.path(rel_out_dir, paste0(prefix, "_05_cox_HR_preservation.csv"))
    )

    p_cox <- ggplot(cox_compare, aes(real_logHR, virtual_logHR)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE) +
      theme_bw() +
      labs(
        title = paste0("Cox logHR preservation: ", set_name),
        subtitle = paste0("r = ", round(cox_cor, 3)),
        x = "Real miRNA Cox logHR",
        y = "Virtual miRNA Cox logHR"
      )

    ggsave(
      file.path(rel_fig_dir, paste0(prefix, "_05_cox_logHR_preservation.pdf")),
      p_cox,
      width = 5,
      height = 4
    )
  }

  # 19.6 subtype preservation ---------------------------------------------
  cluster_res <- list()

  max_k <- min(6, nrow(real_z_sub) - 1)

  if (max_k >= 2 && ncol(real_z_sub) >= 2) {
    for (k in 2:max_k) {
      set.seed(123 + k)

      real_cluster <- tryCatch(
        kmeans(real_z_sub, centers = k, nstart = 50)$cluster,
        error = function(e) NULL
      )

      set.seed(456 + k)

      virtual_cluster <- tryCatch(
        kmeans(virtual_z_sub, centers = k, nstart = 50)$cluster,
        error = function(e) NULL
      )

      if (!is.null(real_cluster) && !is.null(virtual_cluster)) {
        cluster_res[[as.character(k)]] <- tibble(
          miRNA_set = set_name,
          n_miRNA = length(keep_miRNA),
          k = k,
          ARI = ari_score(real_cluster, virtual_cluster),
          NMI = nmi_score(real_cluster, virtual_cluster)
        )
      }
    }
  }

  cluster_df <- bind_rows(cluster_res)

  if (nrow(cluster_df) > 0) {
    write_csv(
      cluster_df,
      file.path(rel_out_dir, paste0(prefix, "_06_subtype_preservation_ARI_NMI.csv"))
    )

    p_ari <- ggplot(cluster_df, aes(k, ARI)) +
      geom_line() +
      geom_point() +
      theme_bw() +
      labs(
        title = paste0("Subtype preservation by ARI: ", set_name),
        x = "Number of clusters",
        y = "ARI"
      )

    ggsave(
      file.path(rel_fig_dir, paste0(prefix, "_06_subtype_preservation_ARI.pdf")),
      p_ari,
      width = 5,
      height = 4
    )

    p_nmi <- ggplot(cluster_df, aes(k, NMI)) +
      geom_line() +
      geom_point() +
      theme_bw() +
      labs(
        title = paste0("Subtype preservation by NMI: ", set_name),
        x = "Number of clusters",
        y = "NMI"
      )

    ggsave(
      file.path(rel_fig_dir, paste0(prefix, "_06_subtype_preservation_NMI.pdf")),
      p_nmi,
      width = 5,
      height = 4
    )
  }

  summary_df <- tibble(
    miRNA_set = set_name,
    n_miRNA = length(keep_miRNA),
    median_sample_correlation = median(sample_cor, na.rm = TRUE),
    mean_sample_correlation = mean(sample_cor, na.rm = TRUE),
    correlation_matrix_preservation = cor_matrix_preservation,
    differential_logFC_correlation = diff_logfc_cor,
    differential_top50_overlap = top_overlap,
    cox_logHR_correlation = cox_cor
  )

  list(
    summary = summary_df,
    subtype = cluster_df,
    sample_correlation = sample_cor_df
  )
}

# 运行四套 downstream 分析 --------------------------------------------------

downstream_results <- list()

for (set_name in names(miRNA_sets)) {
  downstream_results[[set_name]] <- run_downstream_for_miRNA_set(
    set_name = set_name,
    keep_miRNA = miRNA_sets[[set_name]],
    real_z = real_z,
    virtual_z = virtual_z,
    real_expr = real_expr,
    virtual_expr = virtual_expr,
    clinical_test = clinical_test,
    best_model = best_model,
    rel_out_dir = rel_out_dir,
    rel_fig_dir = rel_fig_dir
  )
}

downstream_summary <- bind_rows(lapply(downstream_results, `[[`, "summary")) %>%
  mutate(
    miRNA_set = factor(
      miRNA_set,
      levels = c("All", "Rel_ge_0_5", "Rel_ge_0_6", "Rel_ge_0_7")
    )
  )

write_csv(
  downstream_summary,
  file.path(rel_out_dir, "07_downstream_preservation_summary_by_reliability_threshold.csv")
)

subtype_summary <- bind_rows(lapply(downstream_results, `[[`, "subtype"))

if (nrow(subtype_summary) > 0) {
  subtype_summary <- subtype_summary %>%
    mutate(
      miRNA_set = factor(
        miRNA_set,
        levels = c("All", "Rel_ge_0_5", "Rel_ge_0_6", "Rel_ge_0_7")
      )
    )

  write_csv(
    subtype_summary,
    file.path(rel_out_dir, "08_subtype_preservation_summary_by_reliability_threshold.csv")
  )
}

sample_correlation_all <- bind_rows(lapply(downstream_results, `[[`, "sample_correlation"))

if (nrow(sample_correlation_all) > 0) {
  sample_correlation_all <- sample_correlation_all %>%
    mutate(
      miRNA_set = factor(
        miRNA_set,
        levels = c("All", "Rel_ge_0_5", "Rel_ge_0_6", "Rel_ge_0_7")
      )
    )

  write_csv(
    sample_correlation_all,
    file.path(rel_out_dir, "09_sample_correlation_by_reliability_threshold.csv")
  )

  p_sample_box <- ggplot(
    sample_correlation_all,
    aes(miRNA_set, real_virtual_sample_correlation)
  ) +
    geom_boxplot(outlier.alpha = 0.4) +
    theme_bw() +
    labs(
      title = "Sample-wise real-virtual correlation by reliability threshold",
      x = "miRNA reliability threshold",
      y = "Sample-wise correlation"
    )

  ggsave(
    file.path(rel_fig_dir, "09_sample_correlation_by_reliability_threshold.pdf"),
    p_sample_box,
    width = 6,
    height = 4
  )
}

# 阈值汇总图 ----------------------------------------------------------------

plot_metric_by_threshold <- function(df, metric, ylab, file_name) {
  if (!metric %in% colnames(df)) return(invisible(NULL))
  if (all(is.na(df[[metric]]))) return(invisible(NULL))

  p <- ggplot(df, aes(miRNA_set, .data[[metric]], group = 1)) +
    geom_line() +
    geom_point(size = 2.5) +
    theme_bw() +
    labs(
      title = ylab,
      x = "miRNA reliability threshold",
      y = ylab
    )

  ggsave(
    file.path(rel_fig_dir, file_name),
    p,
    width = 6,
    height = 4
  )
}

plot_metric_by_threshold(
  downstream_summary,
  "median_sample_correlation",
  "Median sample-wise real-virtual correlation",
  "10_summary_median_sample_correlation_by_threshold.pdf"
)

plot_metric_by_threshold(
  downstream_summary,
  "correlation_matrix_preservation",
  "miRNA-miRNA correlation preservation",
  "11_summary_correlation_matrix_preservation_by_threshold.pdf"
)

plot_metric_by_threshold(
  downstream_summary,
  "differential_logFC_correlation",
  "Differential logFC preservation",
  "12_summary_differential_logFC_preservation_by_threshold.pdf"
)

plot_metric_by_threshold(
  downstream_summary,
  "cox_logHR_correlation",
  "Cox logHR preservation",
  "13_summary_cox_logHR_preservation_by_threshold.pdf"
)

if (nrow(subtype_summary) > 0) {
  p_ari_all <- ggplot(subtype_summary, aes(k, ARI, group = miRNA_set, linetype = miRNA_set)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    labs(
      title = "Subtype preservation by ARI across reliability thresholds",
      x = "Number of clusters",
      y = "ARI",
      linetype = "miRNA set"
    )

  ggsave(
    file.path(rel_fig_dir, "14_subtype_ARI_by_reliability_threshold.pdf"),
    p_ari_all,
    width = 6,
    height = 4
  )

  p_nmi_all <- ggplot(subtype_summary, aes(k, NMI, group = miRNA_set, linetype = miRNA_set)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    labs(
      title = "Subtype preservation by NMI across reliability thresholds",
      x = "Number of clusters",
      y = "NMI",
      linetype = "miRNA set"
    )

  ggsave(
    file.path(rel_fig_dir, "15_subtype_NMI_by_reliability_threshold.pdf"),
    p_nmi_all,
    width = 6,
    height = 4
  )
}

# 保存 downstream 单独复跑所需的轻量输入 ----------------------------------

downstream_input <- list(
  best_model = best_model,
  Y_test = Y_test,
  Y_test_raw = Y_test_raw,
  best_pred_test_z = best_pred_test_z,
  best_pred_test_expr = best_pred_test_expr,
  clinical_test = clinical_test,
  miRNA_reliability = miRNA_reliability,
  sample_reliability_test = sample_reliability_test,
  miRNA_sets = miRNA_sets,
  downstream_summary = downstream_summary,
  subtype_summary = subtype_summary
)

saveRDS(
  downstream_input,
  file.path(out_dir, "29_downstream_input_best_model_by_reliability.rds")
)


# 20. 保存模型模块、预测矩阵和 Reli-miR 对象 --------------------------------
# 说明：
#   1. 每个模型单独保存到 output/models。
#   2. 每个模型的 validation/test 预测矩阵单独保存到 output/predictions。
#   3. 最终 Reli-miR 对象只保存轻量 metadata 和路径，不再把所有大模型塞进一个 RDS。
#   4. 后续 downstream 或私有数据预测时，可以按需 readRDS 导入对应模块。

save_relimir_modules <- function(models, pred_val_list, pred_test_list, model_dir, pred_dir) {

  model_registry <- list()
  prediction_registry <- list()

  for (m in names(models)) {
    model_file <- file.path(model_dir, paste0(m, "_model.rds"))
    saveRDS(models[[m]], model_file)
    model_registry[[m]] <- tibble(
      model = m,
      model_file = model_file
    )
  }

  all_pred_models <- union(names(pred_val_list), names(pred_test_list))

  for (m in all_pred_models) {
    val_file <- NA_character_
    test_file <- NA_character_

    if (m %in% names(pred_val_list)) {
      val_file <- file.path(pred_dir, paste0(m, "_prediction_validation_z.rds"))
      saveRDS(pred_val_list[[m]], val_file)
    }

    if (m %in% names(pred_test_list)) {
      test_file <- file.path(pred_dir, paste0(m, "_prediction_test_z.rds"))
      saveRDS(pred_test_list[[m]], test_file)
    }

    prediction_registry[[m]] <- tibble(
      model = m,
      validation_prediction_file = val_file,
      test_prediction_file = test_file
    )
  }

  list(
    model_registry = bind_rows(model_registry),
    prediction_registry = bind_rows(prediction_registry)
  )
}

module_registry <- save_relimir_modules(
  models = models,
  pred_val_list = pred_val_list,
  pred_test_list = pred_test_list,
  model_dir = model_dir,
  pred_dir = pred_dir
)

write_csv(
  module_registry$model_registry,
  file.path(out_dir, "30_model_module_registry.csv")
)

write_csv(
  module_registry$prediction_registry,
  file.path(out_dir, "31_prediction_module_registry.csv")
)

make_model_metadata <- function(models, model_registry) {
  metadata <- list()

  for (m in names(models)) {
    x <- models[[m]]

    # 大对象只保存在独立模型文件中；最终对象里只留轻量信息和路径。
    x$fit <- NULL
    x$pca <- NULL
    x$coef <- NULL

    model_file <- model_registry %>%
      filter(model == m) %>%
      pull(model_file)

    if (length(model_file) == 1) {
      x$model_file <- model_file
    } else {
      x$model_file <- NA_character_
    }

    metadata[[m]] <- x
  }

  metadata
}

model_metadata <- make_model_metadata(
  models = models,
  model_registry = module_registry$model_registry
)

relimir_object <- list(
  project = "Reli-miR",
  version = "final_reliability_thresholds_memory_saving",

  best_model = best_model,
  models = model_metadata,
  model_registry = module_registry$model_registry,
  prediction_registry = module_registry$prediction_registry,

  top_genes = top_genes,
  all_mRNA_genes = colnames(X_train_raw_all),
  miRNA_names = miRNA_names,

  x_median = x_median,
  y_median = y_median,

  x_center = x_center,
  x_scale = x_scale,
  y_center = y_center,
  y_scale = y_scale,

  pathway_object = pathway_object,

  pca_distance = pca_distance,
  n_pc_distance = n_pc_distance,
  distance_center = distance_center,
  distance_variance = distance_variance,
  distance_q95 = distance_q95,

  disagreement_min = disagreement_min,
  disagreement_max = disagreement_max,
  available_prediction_models = available_prediction_models,

  miRNA_level_reliability = miRNA_reliability,
  sample_reliability_validation = sample_reliability_val,
  sample_reliability_test = sample_reliability_test,

  miRNA_sets = miRNA_sets,
  downstream_summary_by_reliability = downstream_summary,
  subtype_summary_by_reliability = subtype_summary,

  global_validation_metric = global_val,
  global_test_metric = global_test,

  split_table = split_table
)

saveRDS(
  relimir_object,
  file.path(out_dir, "32_Reli_miR_model_object_final_lightweight.rds")
)

# 如果后续确实需要所有模型合并在一个对象里，也可以手动保存；默认不保存，避免 RDS 过大。
# saveRDS(models, file.path(out_dir, "32_all_models_heavy_backup.rds"))


# 21. 私有 mRNA 数据预测函数 -------------------------------------------------
# 使用方法：
#   source("Reli_miR_final_reliability_thresholds.R")
#   obj <- readRDS("output/32_Reli_miR_model_object_final_lightweight.rds")
#   private_result <- predict_relimir_private_mrna(private_mrna, obj)
#
# private_mrna 格式：
#   patient_id + mRNA expression matrix
#
# 输出：
#   virtual_miRNA_z
#   virtual_miRNA_expression
#   sample_level_reliability
#   miRNA_level_reliability
#
# 注意：
#   如果 gene_coverage 很低，通常说明 gene symbol / Ensembl ID 不一致。
#   本函数会按需从 obj$model_registry 读取模型文件，不会一次性加载所有模型。

load_relimir_model <- function(obj, model_name) {
  if (!model_name %in% names(obj$models)) {
    stop("对象中不存在模型：", model_name)
  }

  model_file <- obj$models[[model_name]]$model_file

  if (!is.null(model_file) && length(model_file) == 1 && !is.na(model_file) && file.exists(model_file)) {
    return(readRDS(model_file))
  }

  if (!is.null(obj$model_registry)) {
    model_file2 <- obj$model_registry %>%
      filter(model == model_name) %>%
      pull(model_file)

    if (length(model_file2) == 1 && file.exists(model_file2)) {
      return(readRDS(model_file2))
    }
  }

  # 兼容旧版本：如果 obj$models 里仍然直接保存了模型 fit，则直接返回。
  if (!is.null(obj$models[[model_name]]$fit) ||
      !is.null(obj$models[[model_name]]$coef) ||
      !is.null(obj$models[[model_name]]$pca)) {
    return(obj$models[[model_name]])
  }

  stop("找不到模型文件：", model_name, "。请检查 output/models 文件夹是否还在原位置。")
}

predict_relimir_private_mrna <- function(private_mrna, obj, model_name = NULL) {

  if (is.null(model_name)) {
    model_name <- obj$best_model
  }

  if (!"patient_id" %in% colnames(private_mrna)) {
    stop("private_mrna 必须包含 patient_id 列")
  }

  private_mrna <- as.data.frame(private_mrna)

  ids <- as.character(private_mrna$patient_id)

  private_expr <- private_mrna[, setdiff(colnames(private_mrna), "patient_id"), drop = FALSE]
  colnames(private_expr) <- make.unique(colnames(private_expr))

  private_expr[] <- lapply(private_expr, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })

  private_mat0 <- as.matrix(private_expr)
  rownames(private_mat0) <- ids
  storage.mode(private_mat0) <- "numeric"

  # 构建 top gene 矩阵
  required_genes <- obj$top_genes

  common_gene <- intersect(required_genes, colnames(private_mat0))
  gene_coverage <- length(common_gene) / length(required_genes)

  message("Top mRNA gene coverage: ", round(gene_coverage * 100, 2), "%")

  if (gene_coverage < 0.7) {
    warning(
      "Top gene 覆盖率低于 70%，预测结果可能不可靠。请检查 gene symbol / Ensembl ID 是否一致。"
    )
  }

  X_private_raw <- matrix(
    NA_real_,
    nrow = nrow(private_mat0),
    ncol = length(required_genes),
    dimnames = list(ids, required_genes)
  )

  if (length(common_gene) > 0) {
    X_private_raw[, common_gene] <- private_mat0[, common_gene, drop = FALSE]
  }

  med <- obj$x_median[required_genes]
  med[is.na(med)] <- 0

  for (j in seq_len(ncol(X_private_raw))) {
    idx <- !is.finite(X_private_raw[, j]) | is.na(X_private_raw[, j])
    if (any(idx)) {
      X_private_raw[idx, j] <- med[j]
    }
  }

  X_private_z <- sweep(X_private_raw, 2, obj$x_center[required_genes], "-")
  X_private_z <- sweep(X_private_z, 2, obj$x_scale[required_genes], "/")
  X_private_z[!is.finite(X_private_z)] <- 0

  # 如果存在 pathway 模型，构建全部 mRNA 矩阵
  X_private_all_z <- NULL
  all_gene_coverage <- NA_real_

  if (!is.null(obj$pathway_object)) {

    all_required_genes <- obj$all_mRNA_genes
    common_gene_all <- intersect(all_required_genes, colnames(private_mat0))
    all_gene_coverage <- length(common_gene_all) / length(all_required_genes)

    message("All mRNA gene coverage: ", round(all_gene_coverage * 100, 2), "%")

    X_private_all_raw <- matrix(
      NA_real_,
      nrow = nrow(private_mat0),
      ncol = length(all_required_genes),
      dimnames = list(ids, all_required_genes)
    )

    if (length(common_gene_all) > 0) {
      X_private_all_raw[, common_gene_all] <- private_mat0[, common_gene_all, drop = FALSE]
    }

    med_all <- obj$x_median[all_required_genes]
    med_all[is.na(med_all)] <- 0

    for (j in seq_len(ncol(X_private_all_raw))) {
      idx <- !is.finite(X_private_all_raw[, j]) | is.na(X_private_all_raw[, j])
      if (any(idx)) {
        X_private_all_raw[idx, j] <- med_all[j]
      }
    }

    X_private_all_z <- sweep(
      X_private_all_raw,
      2,
      obj$pathway_object$x_all_center[all_required_genes],
      "-"
    )

    X_private_all_z <- sweep(
      X_private_all_z,
      2,
      obj$pathway_object$x_all_scale[all_required_genes],
      "/"
    )

    X_private_all_z[!is.finite(X_private_all_z)] <- 0
  }

  # 对所有可用模型预测，方便计算 agreement
  pred_all_models <- list()

  for (m in obj$available_prediction_models) {

    if (!m %in% names(obj$models)) next

    model <- load_relimir_model(obj, m)

    if (model$type %in% c("ridge", "elastic_net")) {

      pred <- predict(model$fit, newx = X_private_z, s = "lambda.min")
      pred <- to_prediction_matrix(pred, rownames(X_private_z), obj$miRNA_names)

      pred_all_models[[m]] <- pred
    }

    if (model$type == "pls") {

      pred <- X_private_z %*% model$coef
      colnames(pred) <- obj$miRNA_names
      rownames(pred) <- rownames(X_private_z)

      pred_all_models[[m]] <- pred
    }

    if (model$type == "pca_ridge") {

      Z <- predict(model$pca, newdata = X_private_z)[, 1:model$best_k, drop = FALSE]

      pred <- predict(model$fit, newx = Z, s = "lambda.min")
      pred <- to_prediction_matrix(pred, rownames(X_private_z), obj$miRNA_names)

      pred_all_models[[m]] <- pred
    }

    if (model$type == "pca_elasticnet") {

      Z <- predict(model$pca, newdata = X_private_z)[, 1:model$best_k, drop = FALSE]

      pred <- predict(model$fit, newx = Z, s = "lambda.min")
      pred <- to_prediction_matrix(pred, rownames(X_private_z), obj$miRNA_names)

      pred_all_models[[m]] <- pred
    }

    if (model$type == "mlp") {

      Z <- predict(model$pca, newdata = X_private_z)[, 1:model$npc, drop = FALSE]

      pred <- predict(model$fit, Z)
      pred <- as.matrix(pred)
      colnames(pred) <- obj$miRNA_names
      rownames(pred) <- rownames(X_private_z)

      pred_all_models[[m]] <- pred
    }

    if (model$type == "pathway_ridge") {

      if (is.null(obj$pathway_object) || is.null(X_private_all_z)) {
        next
      }

      gene_sets_used <- obj$pathway_object$gene_sets_used

      P_private_raw <- matrix(
        NA_real_,
        nrow = nrow(X_private_all_z),
        ncol = length(gene_sets_used),
        dimnames = list(
          rownames(X_private_all_z),
          names(gene_sets_used)
        )
      )

      for (i in seq_along(gene_sets_used)) {
        gs <- intersect(gene_sets_used[[i]], colnames(X_private_all_z))
        P_private_raw[, i] <- rowMeans(X_private_all_z[, gs, drop = FALSE], na.rm = TRUE)
      }

      P_private_raw[!is.finite(P_private_raw)] <- 0

      P_private_raw <- P_private_raw[, obj$pathway_object$pathway_names, drop = FALSE]

      P_private_z <- sweep(
        P_private_raw,
        2,
        obj$pathway_object$pathway_center[obj$pathway_object$pathway_names],
        "-"
      )

      P_private_z <- sweep(
        P_private_z,
        2,
        obj$pathway_object$pathway_scale[obj$pathway_object$pathway_names],
        "/"
      )

      P_private_z[!is.finite(P_private_z)] <- 0

      pred <- predict(model$fit, newx = P_private_z, s = "lambda.min")
      pred <- to_prediction_matrix(pred, rownames(X_private_z), obj$miRNA_names)

      pred_all_models[[m]] <- pred
    }

    rm(model)
    gc()
  }

  if (!model_name %in% names(pred_all_models)) {
    stop("指定模型没有成功预测：", model_name)
  }

  pred_z <- pred_all_models[[model_name]]

  pred_expr <- sweep(pred_z, 2, obj$y_scale[obj$miRNA_names], "*")
  pred_expr <- sweep(pred_expr, 2, obj$y_center[obj$miRNA_names], "+")

  # sample reliability: distance score
  private_pc <- predict(obj$pca_distance, newdata = X_private_z)[, 1:obj$n_pc_distance, drop = FALSE]

  d_private <- rowSums(
    sweep(private_pc, 2, obj$distance_center, "-")^2 /
      matrix(obj$distance_variance, nrow(private_pc), obj$n_pc_distance, byrow = TRUE)
  )

  distance_score <- exp(-d_private / obj$distance_q95)

  # sample reliability: agreement score
  if (length(pred_all_models) >= 2 &&
      is.finite(obj$disagreement_min) &&
      is.finite(obj$disagreement_max) &&
      abs(obj$disagreement_max - obj$disagreement_min) > 1e-12) {

    pred_array <- array(
      NA_real_,
      dim = c(nrow(pred_z), ncol(pred_z), length(pred_all_models)),
      dimnames = list(
        rownames(pred_z),
        colnames(pred_z),
        names(pred_all_models)
      )
    )

    for (i in seq_along(pred_all_models)) {
      pred_array[, , i] <- pred_all_models[[i]]
    }

    model_disagreement <- apply(pred_array, 1, function(x) {
      mean(apply(x, 1, sd, na.rm = TRUE), na.rm = TRUE)
    })

    agreement_score <- 1 - (model_disagreement - obj$disagreement_min) /
      (obj$disagreement_max - obj$disagreement_min)

    agreement_score <- pmin(pmax(agreement_score, 0), 1)

  } else {

    model_disagreement <- rep(NA_real_, nrow(pred_z))
    agreement_score <- rep(NA_real_, nrow(pred_z))
  }

  sample_level_reliability <- rowMeans(
    cbind(distance_score, agreement_score),
    na.rm = TRUE
  )

  sample_reliability <- tibble(
    patient_id = rownames(pred_z),
    gene_coverage = gene_coverage,
    all_gene_coverage = all_gene_coverage,
    distance_score = as.numeric(distance_score),
    model_disagreement = as.numeric(model_disagreement),
    agreement_score = as.numeric(agreement_score),
    sample_level_reliability = as.numeric(sample_level_reliability)
  )

  list(
    model = model_name,
    gene_coverage = gene_coverage,
    all_gene_coverage = all_gene_coverage,
    virtual_miRNA_z = pred_z,
    virtual_miRNA_expression = pred_expr,
    sample_level_reliability = sample_reliability,
    miRNA_level_reliability = obj$miRNA_level_reliability,
    all_model_predictions_z = pred_all_models
  )
}


# 22. 保存 sessionInfo -------------------------------------------------------

sink(file.path(out_dir, "33_sessionInfo.txt"))
print(sessionInfo())
sink()


# 23. 最终总结 --------------------------------------------------------------

best_test_row <- global_test %>%
  filter(model == best_model)

final_summary <- tibble(
  item = c(
    "n_train",
    "n_validation",
    "n_test",
    "n_selected_mRNA",
    "n_miRNA",
    "best_model_selected_by_validation",
    "best_model_test_mean_R2",
    "best_model_test_median_Pearson",
    "best_model_test_mean_RMSE",
    "miRNA_reliability_vs_test_error_spearman",
    "sample_reliability_vs_test_error_spearman",
    "n_miRNA_All",
    "n_miRNA_Rel_ge_0_5",
    "n_miRNA_Rel_ge_0_6",
    "n_miRNA_Rel_ge_0_7"
  ),
  value = c(
    nrow(X_train),
    nrow(X_val),
    nrow(X_test),
    ncol(X_train),
    ncol(Y_train),
    best_model,
    best_test_row$mean_R2,
    best_test_row$median_Pearson,
    best_test_row$mean_RMSE,
    miRNA_rel_test_error_cor,
    sample_rel_test_error_cor,
    length(miRNA_sets$All),
    length(miRNA_sets$Rel_ge_0_5),
    length(miRNA_sets$Rel_ge_0_6),
    length(miRNA_sets$Rel_ge_0_7)
  )
)

write_csv(
  final_summary,
  file.path(out_dir, "34_Reli_miR_final_summary.csv")
)

message("\nReli-miR final 分析完成。")
message("最佳模型：", best_model)
message("结果保存位置：", out_dir)
message("模型模块保存位置：", model_dir)
message("预测矩阵保存位置：", pred_dir)
message("Reliability 阈值 downstream 结果保存位置：", rel_out_dir)
