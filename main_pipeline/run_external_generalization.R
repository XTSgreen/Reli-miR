options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(pls)
  library(nnet)
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

normalize_expression_df <- function(df, id_col = "patient_id") {
  df <- as.data.frame(df)
  stopifnot(id_col %in% colnames(df))

  df %>%
    mutate(patient_id = as.character(.data[[id_col]])) %>%
    distinct(patient_id, .keep_all = TRUE) %>%
    { .[, c(id_col, setdiff(colnames(.), id_col)), drop = FALSE] }
}

build_external_matrices <- function(external_mrna, external_mirna, obj, id_col = "patient_id") {
  external_mrna <- normalize_expression_df(external_mrna, id_col = id_col)
  external_mirna <- normalize_expression_df(external_mirna, id_col = id_col)

  common_id <- intersect(external_mrna[[id_col]], external_mirna[[id_col]])
  if (length(common_id) == 0) {
    stop("外部 mRNA 与 miRNA 没有可对齐的 patient_id。")
  }

  mrna2 <- external_mrna[match(common_id, external_mrna[[id_col]]), , drop = FALSE]
  mirna2 <- external_mirna[match(common_id, external_mirna[[id_col]]), , drop = FALSE]

  x_expr <- mrna2[, setdiff(colnames(mrna2), id_col), drop = FALSE]
  y_expr <- mirna2[, setdiff(colnames(mirna2), id_col), drop = FALSE]
  colnames(x_expr) <- make.unique(colnames(x_expr))
  colnames(y_expr) <- make.unique(colnames(y_expr))

  x_expr[] <- lapply(x_expr, function(x) suppressWarnings(as.numeric(as.character(x))))
  y_expr[] <- lapply(y_expr, function(x) suppressWarnings(as.numeric(as.character(x))))

  x_mat0 <- as.matrix(x_expr)
  y_mat0 <- as.matrix(y_expr)
  rownames(x_mat0) <- common_id
  rownames(y_mat0) <- common_id
  storage.mode(x_mat0) <- "numeric"
  storage.mode(y_mat0) <- "numeric"

  required_genes <- obj$top_genes
  common_gene <- intersect(required_genes, colnames(x_mat0))
  gene_coverage <- length(common_gene) / length(required_genes)

  x_private_raw <- matrix(
    NA_real_,
    nrow = nrow(x_mat0),
    ncol = length(required_genes),
    dimnames = list(common_id, required_genes)
  )
  if (length(common_gene) > 0) {
    x_private_raw[, common_gene] <- x_mat0[, common_gene, drop = FALSE]
  }

  med <- obj$x_median[required_genes]
  med[is.na(med)] <- 0
  for (j in seq_len(ncol(x_private_raw))) {
    idx <- !is.finite(x_private_raw[, j]) | is.na(x_private_raw[, j])
    if (any(idx)) {
      x_private_raw[idx, j] <- med[j]
    }
  }

  x_private_z <- sweep(x_private_raw, 2, obj$x_center[required_genes], "-")
  x_private_z <- sweep(x_private_z, 2, obj$x_scale[required_genes], "/")
  x_private_z[!is.finite(x_private_z)] <- 0

  x_private_all_z <- NULL
  all_gene_coverage <- NA_real_
  if (!is.null(obj$all_mRNA_genes) && !is.null(obj$pathway_object)) {
    all_required_genes <- obj$all_mRNA_genes
    common_gene_all <- intersect(all_required_genes, colnames(x_mat0))
    all_gene_coverage <- length(common_gene_all) / length(all_required_genes)

    x_private_all_raw <- matrix(
      NA_real_,
      nrow = nrow(x_mat0),
      ncol = length(all_required_genes),
      dimnames = list(common_id, all_required_genes)
    )
    if (length(common_gene_all) > 0) {
      x_private_all_raw[, common_gene_all] <- x_mat0[, common_gene_all, drop = FALSE]
    }

    med_all <- obj$x_median[all_required_genes]
    med_all[is.na(med_all)] <- 0
    for (j in seq_len(ncol(x_private_all_raw))) {
      idx <- !is.finite(x_private_all_raw[, j]) | is.na(x_private_all_raw[, j])
      if (any(idx)) {
        x_private_all_raw[idx, j] <- med_all[j]
      }
    }

    center_all <- obj$pathway_object$x_all_center[all_required_genes]
    scale_all <- obj$pathway_object$x_all_scale[all_required_genes]
    x_private_all_z <- sweep(x_private_all_raw, 2, center_all, "-")
    x_private_all_z <- sweep(x_private_all_z, 2, scale_all, "/")
    x_private_all_z[!is.finite(x_private_all_z)] <- 0
  }

  common_mirna <- intersect(obj$miRNA_names, colnames(y_mat0))
  if (length(common_mirna) == 0) {
    stop("外部 miRNA 与训练对象没有共同 miRNA。")
  }

  y_true_z <- y_mat0[, common_mirna, drop = FALSE]
  y_true_z <- sweep(y_true_z, 2, obj$y_center[common_mirna], "-")
  y_true_z <- sweep(y_true_z, 2, obj$y_scale[common_mirna], "/")
  y_true_z[!is.finite(y_true_z)] <- NA_real_

  list(
    patient_ids = common_id,
    x_private_z = x_private_z,
    x_private_all_z = x_private_all_z,
    y_true_z = y_true_z,
    common_mirna = common_mirna,
    gene_coverage = gene_coverage,
    all_gene_coverage = all_gene_coverage
  )
}

predict_external_with_model <- function(model, obj, x_private_z, x_private_all_z = NULL) {
  if (model$type %in% c("ridge", "elastic_net")) {
    pred <- predict(model$fit, newx = x_private_z, s = "lambda.min")
    return(to_prediction_matrix(pred, rownames(x_private_z), obj$miRNA_names))
  }

  if (model$type == "pls") {
    pred <- x_private_z %*% model$coef
    pred <- as.matrix(pred)
    rownames(pred) <- rownames(x_private_z)
    colnames(pred) <- obj$miRNA_names
    return(pred)
  }

  if (model$type == "pca_ridge") {
    z <- predict(model$pca, newdata = x_private_z)[, 1:model$best_k, drop = FALSE]
    pred <- predict(model$fit, newx = z, s = "lambda.min")
    return(to_prediction_matrix(pred, rownames(x_private_z), obj$miRNA_names))
  }

  if (model$type == "pca_elasticnet") {
    z <- predict(model$pca, newdata = x_private_z)[, 1:model$best_k, drop = FALSE]
    pred <- predict(model$fit, newx = z, s = "lambda.min")
    return(to_prediction_matrix(pred, rownames(x_private_z), obj$miRNA_names))
  }

  if (model$type == "mlp") {
    z <- predict(model$pca, newdata = x_private_z)[, 1:model$npc, drop = FALSE]
    pred <- predict(model$fit, z)
    pred <- as.matrix(pred)
    rownames(pred) <- rownames(x_private_z)
    colnames(pred) <- obj$miRNA_names
    return(pred)
  }

  if (model$type == "pathway_ridge") {
    if (is.null(obj$pathway_object) || is.null(x_private_all_z)) {
      stop("Pathway_Ridge 需要 obj$pathway_object 和全基因矩阵。")
    }

    gene_sets_used <- obj$pathway_object$gene_sets_used
    p_private_raw <- matrix(
      NA_real_,
      nrow = nrow(x_private_all_z),
      ncol = length(gene_sets_used),
      dimnames = list(rownames(x_private_all_z), names(gene_sets_used))
    )

    for (i in seq_along(gene_sets_used)) {
      gs <- intersect(gene_sets_used[[i]], colnames(x_private_all_z))
      p_private_raw[, i] <- rowMeans(x_private_all_z[, gs, drop = FALSE], na.rm = TRUE)
    }

    p_private_raw[!is.finite(p_private_raw)] <- 0
    p_private_raw <- p_private_raw[, obj$pathway_object$pathway_names, drop = FALSE]

    p_private_z <- sweep(
      p_private_raw,
      2,
      obj$pathway_object$pathway_center[obj$pathway_object$pathway_names],
      "-"
    )
    p_private_z <- sweep(
      p_private_z,
      2,
      obj$pathway_object$pathway_scale[obj$pathway_object$pathway_names],
      "/"
    )
    p_private_z[!is.finite(p_private_z)] <- 0

    pred <- predict(model$fit, newx = p_private_z, s = "lambda.min")
    return(to_prediction_matrix(pred, rownames(x_private_z), obj$miRNA_names))
  }

  stop("不支持的模型类型：", model$type)
}

evaluate_external_models <- function(external_mrna,
                                     external_mirna,
                                     obj,
                                     model_files,
                                     dataset_name = "external") {
  built <- build_external_matrices(external_mrna, external_mirna, obj)

  metric_list <- list()
  global_list <- list()
  coverage_list <- list()
  prediction_list <- list()

  for (model_name in names(model_files)) {
    model_file <- model_files[[model_name]]
    if (!file.exists(model_file)) next

    model <- readRDS(model_file)
    pred_z <- tryCatch(
      predict_external_with_model(
        model = model,
        obj = obj,
        x_private_z = built$x_private_z,
        x_private_all_z = built$x_private_all_z
      ),
      error = function(e) {
        warning("模型 ", model_name, " 外部预测失败：", conditionMessage(e))
        NULL
      }
    )
    if (is.null(pred_z)) next

    pred_common <- pred_z[, built$common_mirna, drop = FALSE]
    metric_df <- build_metric_df(
      real = built$y_true_z,
      pred = pred_common,
      dataset = dataset_name,
      model_name = model_name
    )
    global_df <- build_global_df(metric_df)

    metric_list[[model_name]] <- metric_df
    global_list[[model_name]] <- global_df
    prediction_list[[model_name]] <- pred_common
    coverage_list[[model_name]] <- tibble(
      model = model_name,
      n_samples = nrow(pred_common),
      n_miRNA_overlap = length(built$common_mirna),
      gene_coverage = built$gene_coverage,
      all_gene_coverage = built$all_gene_coverage
    )
  }

  if (length(global_list) == 0) {
    stop("没有任何模型成功完成外部泛化评估。")
  }

  global_metrics <- bind_rows(global_list) %>%
    arrange(desc(median_Pearson), desc(mean_R2), mean_RMSE)

  list(
    global_metrics = global_metrics,
    metrics_by_miRNA = bind_rows(metric_list),
    coverage_summary = bind_rows(coverage_list),
    common_miRNA = built$common_mirna,
    patient_ids = built$patient_ids,
    predictions_z = prediction_list,
    best_model_by_median_pearson = global_metrics$model[[1]]
  )
}

load_external_data_rda <- function(data_rda_path) {
  env <- new.env(parent = emptyenv())
  load(data_rda_path, envir = env)

  if (!exists("mrna", envir = env) || !exists("mirna", envir = env)) {
    stop("外部 data.rda 里必须包含 mrna 和 mirna 对象。")
  }

  list(
    clinical = if (exists("clinical", envir = env)) get("clinical", envir = env) else NULL,
    mrna = get("mrna", envir = env),
    mirna = get("mirna", envir = env)
  )
}

discover_model_files <- function(project_dir) {
  model_dir <- file.path(project_dir, "output", "models")
  files <- list.files(model_dir, pattern = "_model\\.rds$", full.names = TRUE)
  names(files) <- sub("_model\\.rds$", "", basename(files))
  files
}

parse_args <- function(args) {
  parsed <- list(
    project_dir = normalizePath(".", winslash = "/", mustWork = TRUE),
    external_data = "data/livergabage_cancer_data/chol/data.rda",
    object_rds = "output/32_Reli_miR_model_object_final_lightweight.rds",
    output_dir = "output/external_generalization"
  )

  if (length(args) %% 2 != 0) {
    stop("参数格式错误。请使用 --key value。")
  }

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    value <- args[[i + 1]]
    if (key == "--project-dir") parsed$project_dir <- value
    if (key == "--external-data") parsed$external_data <- value
    if (key == "--object-rds") parsed$object_rds <- value
    if (key == "--output-dir") parsed$output_dir <- value
    i <- i + 2
  }

  parsed$project_dir <- normalizePath(parsed$project_dir, winslash = "/", mustWork = TRUE)
  parsed$external_data <- normalizePath(file.path(parsed$project_dir, parsed$external_data), winslash = "/", mustWork = TRUE)
  parsed$object_rds <- normalizePath(file.path(parsed$project_dir, parsed$object_rds), winslash = "/", mustWork = TRUE)
  parsed$output_dir <- file.path(parsed$project_dir, parsed$output_dir)
  parsed
}

run_external_generalization <- function(project_dir,
                                        external_data,
                                        object_rds,
                                        output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  obj <- readRDS(object_rds)
  external <- load_external_data_rda(external_data)
  model_files <- discover_model_files(project_dir)

  result <- evaluate_external_models(
    external_mrna = external$mrna,
    external_mirna = external$mirna,
    obj = obj,
    model_files = model_files,
    dataset_name = basename(dirname(external_data))
  )

  readr::write_csv(
    result$global_metrics,
    file.path(output_dir, "external_global_metrics.csv")
  )
  readr::write_csv(
    result$metrics_by_miRNA,
    file.path(output_dir, "external_metrics_by_miRNA.csv")
  )
  readr::write_csv(
    result$coverage_summary,
    file.path(output_dir, "external_gene_coverage_summary.csv")
  )

  summary_list <- list(
    external_data = external_data,
    best_model_by_external_median_pearson = result$best_model_by_median_pearson,
    n_models = nrow(result$global_metrics),
    n_samples = length(result$patient_ids),
    n_miRNA_overlap = length(result$common_miRNA)
  )

  write_json(
    summary_list,
    file.path(output_dir, "external_generalization_summary.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )

  invisible(result)
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  run_external_generalization(
    project_dir = args$project_dir,
    external_data = args$external_data,
    object_rds = args$object_rds,
    output_dir = args$output_dir
  )
}
