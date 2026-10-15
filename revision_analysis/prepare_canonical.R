#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(jsonlite))
out <- "new/results/E0"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
load_env <- new.env(parent = emptyenv())
load("data/data.rda", envir = load_env)
obj_path <- "output/32_Reli_miR_model_object_final_lightweight.rds"
obj <- readRDS(obj_path)
model_path <- "output/models/PCA_ElasticNet_model.rds"
model <- readRDS(model_path)
read_mat <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  m <- as.matrix(d[, setdiff(names(d), "patient_id"), drop = FALSE])
  rownames(m) <- as.character(d$patient_id)
  storage.mode(m) <- "double"
  m
}
inp <- "output/gpu_model_inputs"
can <- list()
for (part in c("train", "cal", "test")) {
  old <- if (part == "cal") "validation" else part
  can[[paste0("X_", part, "_z")]] <- read_mat(file.path(inp, paste0(old, "_x_z.csv")))
  can[[paste0("Y_", part, "_z")]] <- read_mat(file.path(inp, paste0(old, "_y_z.csv")))
}
can$gene_names <- colnames(can$X_train_z)
can$miRNA_names <- colnames(can$Y_train_z)
can$x_center <- obj$x_center[can$gene_names]
can$x_scale <- obj$x_scale[can$gene_names]
can$y_center <- obj$y_center[can$miRNA_names]
can$y_scale <- obj$y_scale[can$miRNA_names]
stopifnot(!anyNA(c(can$x_center, can$x_scale, can$y_center, can$y_scale)))
can$x_median <- obj$x_median[can$gene_names]
can$y_median <- obj$y_median[can$miRNA_names]
clinical <- load_env$clinical
mrna <- load_env$mrna
mirna <- load_env$mirna
checks <- list()
for (part in c("train", "cal", "test")) {
  x <- can[[paste0("X_", part, "_z")]]
  y <- can[[paste0("Y_", part, "_z")]]
  stopifnot(identical(rownames(x), rownames(y)),
            identical(colnames(x), can$gene_names), identical(colnames(y), can$miRNA_names))
  ids <- rownames(y)
  xx <- as.matrix(mrna[match(ids, mrna$patient_id), can$gene_names, drop = FALSE])
  yy <- as.matrix(mirna[match(ids, mirna$patient_id), can$miRNA_names, drop = FALSE])
  storage.mode(xx) <- "double"; storage.mode(yy) <- "double"
  for (j in seq_len(ncol(xx))) xx[!is.finite(xx[,j]),j] <- can$x_median[j]
  for (j in seq_len(ncol(yy))) yy[!is.finite(yy[,j]),j] <- can$y_median[j]
  rownames(xx) <- rownames(yy) <- ids
  expected_x <- sweep(sweep(xx, 2, can$x_center, "-"), 2, can$x_scale, "/")
  expected_y <- sweep(sweep(yy, 2, can$y_center, "-"), 2, can$y_scale, "/")
  errx <- max(abs(expected_x - x)); erry <- max(abs(expected_y - y))
  stopifnot(errx < 1e-9, erry < 1e-9)
  checks[[length(checks)+1]] <- data.frame(check=paste0(part,"_export_vs_raw_truth"),
    max_abs_x_difference=errx,max_abs_y_difference=erry,passed=TRUE)
  can[[paste0("X_",part,"_native")]] <- xx
  can[[paste0("Y_",part,"_native")]] <- yy
  aligned_clinical <- clinical[match(ids, clinical$patient_id),,drop=FALSE]
  rownames(aligned_clinical) <- ids
  stopifnot(identical(as.character(aligned_clinical$patient_id), ids))
  can[[paste0("clinical_",part)]] <- aligned_clinical
}
stopifnot(nrow(can$X_train_z)==222,nrow(can$X_cal_z)==72,nrow(can$X_test_z)==72)
stopifnot(length(intersect(rownames(can$X_train_z),rownames(can$X_cal_z)))==0,
          length(intersect(rownames(can$X_train_z),rownames(can$X_test_z)))==0,
          length(intersect(rownames(can$X_cal_z),rownames(can$X_test_z)))==0)
as_pred <- function(p, rows, cols) {
  if (length(dim(p))==3) p <- p[,,1,drop=FALSE]
  p <- matrix(p,nrow=length(rows),ncol=length(cols),dimnames=list(rows,cols))
  p
}
metrics <- function(y,p) {
  mse <- colMeans((y-p)^2)
  r2 <- 1-colSums((y-p)^2)/colSums(sweep(y,2,colMeans(y),"-")^2)
  pearson <- vapply(seq_len(ncol(y)),function(j){
    if(sd(y[,j])<=1e-12 || sd(p[,j])<=1e-12) return(NA_real_)
    cor(y[,j],p[,j])
  },numeric(1))
  data.frame(miRNA=colnames(y),R2=r2,Pearson=pearson,RMSE=sqrt(mse))
}
for (part in c("cal","test")) {
  old <- if(part=="cal") "validation" else part
  saved <- readRDS(paste0("output/predictions/PCA_ElasticNet_prediction_",old,"_z.rds"))
  truth <- can[[paste0("Y_",part,"_z")]]
  stopifnot(setequal(rownames(saved),rownames(truth)),setequal(colnames(saved),colnames(truth)))
  saved <- saved[rownames(truth),colnames(truth),drop=FALSE]
  z <- predict(model$pca,newdata=can[[paste0("X_",part,"_z")]])[,seq_len(model$best_k),drop=FALSE]
  current <- as_pred(predict(model$fit,newx=z,s="lambda.min"),rownames(truth),colnames(truth))
  discrepancy <- max(abs(saved-current))
  stopifnot(discrepancy<1e-9)
  can[[paste0("pred_",part,"_z")]] <- current
  met <- metrics(truth,current)
  write.csv(met,file.path(out,paste0("PCA_ElasticNet_",part,"_metrics.csv")),row.names=FALSE)
  cat(part,"mean_R2",mean(met$R2),"median_Pearson",median(met$Pearson,na.rm=TRUE),
      "mean_RMSE",mean(met$RMSE),"saved_vs_model",discrepancy,"\n")
}
can$score_table <- metrics(can$Y_cal_z,can$pred_cal_z)
clamp <- function(x) pmax(0,pmin(1,x))
can$score_table$Rel <- 0.5*clamp(can$score_table$R2)+0.5*clamp(can$score_table$Pearson)
can$score_table$Rel[!is.finite(can$score_table$Rel)] <- 0
write.csv(can$score_table,file.path(out,"PCA_ElasticNet_reliability_from_cal.csv"),row.names=FALSE)
counts <- do.call(rbind,lapply(c(.5,.6,.7),function(t){
  ss <- can$score_table[can$score_table$Rel>=t,,drop=FALSE]
  write.csv(ss,file.path(out,paste0("miRNA_Rel_ge_",t,".csv")),row.names=FALSE)
  data.frame(threshold=t,retained=nrow(ss),total=nrow(can$score_table),coverage=nrow(ss)/nrow(can$score_table))
}))
test_met <- metrics(can$Y_test_z,can$pred_test_z)
rho <- cor(can$score_table$Rel,test_met$RMSE,method="spearman")
can$metadata_model_path <- model_path
can$lightweight_object_path <- obj_path
can$provenance <- list(raw_data="data/data.rda",inputs=inp,model=model_path,
  predictions="recomputed model predictions verified against output/predictions",
  fixed_split="original exploratory split; calibration used in original tuning",
  data_source_status="GDC case projects verified: LIHC 366 and CHOL 35; primary mRNA matches official Xena HiSeqV2 in all 6662664 cells",
  expression_unit="mRNA log2(norm_count+1), RSEM normalized counts; miRNA log2(RPM+1) supported by five GDC case files",
  mRNA_source="https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/HiSeqV2.gz",
  mRNA_metadata_version="2017-10-13", OS_time_unit="days; all 366 values exactly match GDC death or last-follow-up days",
  provenance_evidence="new/results/E0/public_provenance_search.md",seed_original=123)
saveRDS(can,file.path(out,"canonical_input.rds"))
write.csv(do.call(rbind,checks),file.path(out,"raw_truth_alignment_checks.csv"),row.names=FALSE)
write.csv(counts,file.path(out,"recomputed_threshold_counts.csv"),row.names=FALSE)
old17 <- read.csv("output/17_miRNA_level_reliability_from_validation.csv",check.names=FALSE)
old_names <- names(old17)
cmp <- data.frame(source=c("canonical","original_manuscript","stale_summary_csv"),
  model=c("PCA_ElasticNet","PCA_ElasticNet","PCA_Ridge"),
  miRNA_reliability_error_rho=c(rho,-.5978,-.590073123394931))
write.csv(cmp,file.path(out,"version_comparison.csv"),row.names=FALSE)
audit <- list(technical_gate="passed",n_train=222,n_cal=72,n_test=72,n_genes=ncol(can$X_train_z),
  n_miRNA=ncol(can$Y_train_z),model=model$type,k=model$best_k,alpha=model$best_alpha,
  lambda_min=model$fit$lambda.min,reliability_test_error_spearman=rho,
  thresholds=counts,original_score_columns=old_names,
  remaining_metadata=c("original project download date and preprocessing scripts unavailable; full primary mRNA numerical match and units verified; miRNA five-case verification supports transformation but not all-case download versions",
    "original reliability thresholds not proven prospectively fixed"))
write_json(audit,file.path(out,"E0_audit.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
cat("Canonical technical gate PASSED; rho",rho,"threshold counts",paste(counts$retained,collapse=","),"\n")
