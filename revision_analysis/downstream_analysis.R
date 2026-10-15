# E3 downstream analyses for the canonical PCA_ElasticNet input.
# Run only after new/results/E0/canonical_input.rds is frozen and reviewed.
# Usage: Rscript new/code/downstream_analysis.R [canonical_input.rds] [output_dir]

options(stringsAsFactors = FALSE)

ari_score <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  if (n < 2) return(NA_real_)
  c2 <- function(z) z * (z - 1) / 2
  idx <- sum(c2(tab)); a <- sum(c2(rowSums(tab))); b <- sum(c2(colSums(tab)))
  exp_idx <- a * b / c2(n); max_idx <- (a + b) / 2
  if (abs(max_idx - exp_idx) < 1e-12) return(NA_real_)
  (idx - exp_idx) / (max_idx - exp_idx)
}

nmi_score <- function(x, y) {
  tab <- table(x, y); p <- tab / sum(tab); px <- rowSums(p); py <- colSums(p)
  nz <- p > 0
  if (!any(nz)) return(NA_real_)
  rr <- row(p)[nz]; cc <- col(p)[nz]
  mi <- sum(p[nz] * log(p[nz] / (px[rr] * py[cc])))
  hx <- -sum(px[px > 0] * log(px[px > 0])); hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx < 1e-12 || hy < 1e-12) return(NA_real_)
  mi / sqrt(hx * hy)
}

as_numeric_matrix <- function(x, label) {
  x <- as.matrix(x); storage.mode(x) <- "double"
  if (is.null(rownames(x)) || is.null(colnames(x))) stop(label, " must have patient and feature names")
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop(label, " has duplicate names")
  x
}

stage_group <- function(clinical, ids) {
  if (is.null(clinical) || !"stage" %in% names(clinical)) return(NULL)
  if (nrow(clinical) != length(ids)) stop("Clinical rows cannot be aligned to patients")
  raw <- toupper(gsub("\\s+", "", as.character(clinical$stage)))
  raw <- sub("^STAGE", "", raw)
  early <- raw %in% c("I", "IA", "IB", "IC", "II", "IIA", "IIB", "IIC", "1", "1A", "1B", "1C", "2", "2A", "2B", "2C")
  late <- raw %in% c("III", "IIIA", "IIIB", "IIIC", "IV", "IVA", "IVB", "IVC", "3", "3A", "3B", "3C", "4", "4A", "4B", "4C")
  out <- rep(NA_character_, length(raw)); out[early] <- "Early"; out[late] <- "Late"
  if (sum(out == "Early", na.rm = TRUE) >= 5L && sum(out == "Late", na.rm = TRUE) >= 5L) {
    attr(out, "source") <- "stage_Early_vs_Late"
    return(out)
  }
  NULL
}

event_group <- function(clinical, ids) {
  if (is.null(clinical) || !"OS_status" %in% names(clinical)) return(NULL)
  status <- suppressWarnings(as.integer(as.character(clinical$OS_status)))
  if (!all(na.omit(unique(status)) %in% c(0L, 1L))) return(NULL)
  out <- ifelse(status == 1L, "Event", ifelse(status == 0L, "No_event", NA_character_))
  if (sum(out == "Event", na.rm = TRUE) >= 5L && sum(out == "No_event", na.rm = TRUE) >= 5L) {
    attr(out, "source") <- "OS_status_event_vs_no_event"
    return(out)
  }
  NULL
}

align_clinical <- function(clinical, ids, label = "clinical") {
  if (is.null(clinical) || !"patient_id" %in% names(clinical)) stop(label, " must include patient_id for explicit alignment")
  key <- as.character(clinical$patient_id)
  if (anyNA(key) || anyDuplicated(key)) stop(label, "$patient_id must be complete and unique")
  idx <- match(ids, key)
  if (anyNA(idx)) stop(label, " is missing patient_id(s): ", paste(ids[is.na(idx)], collapse = ","))
  out <- clinical[idx, , drop = FALSE]
  rownames(out) <- ids
  out
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

feature_stats <- function(M, group) {
  observed <- if (is.null(group)) character() else unique(group[!is.na(group)])
  lev <- if (all(c("Early","Late") %in% observed)) c("Early","Late") else if (all(c("No_event","Event") %in% observed)) c("No_event","Event") else character()
  if (length(lev) != 2L) return(data.frame(miRNA = colnames(M), effect = NA_real_, p = NA_real_))
  i1 <- which(group == lev[1]); i2 <- which(group == lev[2])
  eff <- p <- rep(NA_real_, ncol(M))
  for (j in seq_len(ncol(M))) {
    a <- M[i1, j]; b <- M[i2, j]
    if (length(a) < 5L || length(b) < 5L || !all(is.finite(c(a, b)))) next
    # Y and predictions are already expressed in the common training-SD scale.
    eff[j] <- mean(b) - mean(a)
    p[j] <- tryCatch(t.test(b, a, var.equal = FALSE)$p.value, error = function(e) NA_real_)
  }
  data.frame(miRNA = colnames(M), effect = eff, p = p, stringsAsFactors = FALSE)
}

cox_stats <- function(M, clinical, ids) {
  empty <- data.frame(miRNA = colnames(M), logHR = NA_real_, p = NA_real_, n = 0L, events = 0L, censored = 0L)
  if (is.null(clinical) || !all(c("OS_time", "OS_status") %in% names(clinical))) return(empty)
  time <- suppressWarnings(as.numeric(as.character(clinical$OS_time)))
  status <- suppressWarnings(as.integer(as.character(clinical$OS_status)))
  for (j in seq_len(ncol(M))) {
    ok <- is.finite(time) & time > 0 & status %in% c(0L, 1L) & is.finite(M[, j])
    n <- sum(ok); ev <- sum(status[ok] == 1L)
    empty$n[j] <- n; empty$events[j] <- ev; empty$censored[j] <- n - ev
    if (n < 10 || ev < 3 || sd(M[ok, j]) == 0) next
    # Y and predictions are already training-SD standardized; never divide by native y_scale again.
    fit <- tryCatch(survival::coxph(survival::Surv(time[ok], status[ok]) ~ M[ok, j], ties = "efron"), error = function(e) NULL)
    if (!is.null(fit)) {
      z <- tryCatch(summary(fit), error = function(e) NULL)
      if (!is.null(z) && nrow(z$coefficients) == 1) {
        empty$logHR[j] <- unname(z$coefficients[1, "coef"])
        empty$p[j] <- unname(z$coefficients[1, "Pr(>|z|)"])
      }
    }
  }
  empty
}

de_overlap <- function(real_de, pred_de, m, k) {
  rf <- is.finite(real_de$effect); pf <- is.finite(pred_de$effect)
  nr <- sum(rf); np <- sum(pf)
  if (m <= k || nr < k || np < k) return(c(overlap=NA_real_,real_n=nr,pred_n=np))
  r <- real_de$miRNA[rf][order(-abs(real_de$effect[rf]), real_de$miRNA[rf])][seq_len(k)]
  p <- pred_de$miRNA[pf][order(-abs(pred_de$effect[pf]), pred_de$miRNA[pf])][seq_len(k)]
  c(overlap=length(intersect(r,p))/k,real_n=nr,pred_n=np)
}

procrustes_disparity <- function(A, B) {
  if (nrow(A) < 3 || ncol(A) < 1 || ncol(B) < 1) return(NA_real_)
  # Center and normalize both patient-score configurations, then optimally rotate B.
  A <- scale(A, center = TRUE, scale = FALSE); B <- scale(B, center = TRUE, scale = FALSE)
  if (sqrt(sum(A^2)) == 0 || sqrt(sum(B^2)) == 0) return(NA_real_)
  A <- A / sqrt(sum(A^2)); B <- B / sqrt(sum(B^2))
  q <- min(ncol(A), ncol(B)); A <- A[, seq_len(q), drop = FALSE]; B <- B[, seq_len(q), drop = FALSE]
  sv <- svd(crossprod(B, A)); rot <- sv$u %*% t(sv$v)
  mean(rowSums((A - B %*% rot)^2))
}

pca_metrics <- function(A, B) {
  # Common projection uses the real-data PCA basis and its center/scale for both matrices.
  k <- min(10L, ncol(A), nrow(A) - 1L)
  if (k < 1L) return(c(common_projection_disparity = NA_real_, procrustes_disparity = NA_real_))
  pc <- tryCatch(prcomp(A, center = TRUE, scale. = FALSE, rank. = k), error = function(e) NULL)
  if (is.null(pc)) return(c(common_projection_disparity = NA_real_, procrustes_disparity = NA_real_))
  Ac <- sweep(A, 2, pc$center, "-") %*% pc$rotation[, seq_len(k), drop = FALSE]
  Bc <- sweep(B, 2, pc$center, "-") %*% pc$rotation[, seq_len(k), drop = FALSE]
  common <- mean(rowSums((scale(Ac, center = TRUE, scale = FALSE) - scale(Bc, center = TRUE, scale = FALSE))^2))
  pb <- tryCatch(prcomp(B, center = TRUE, scale. = FALSE, rank. = k), error = function(e) NULL)
  pro <- if (is.null(pb)) NA_real_ else procrustes_disparity(pc$x[, seq_len(min(k, ncol(pc$x))), drop = FALSE], pb$x[, seq_len(min(k, ncol(pb$x))), drop = FALSE])
  c(common_projection_disparity = common, procrustes_disparity = pro)
}

cluster_metrics <- function(A, B, seed = 101L, nstart = 50L) {
  out <- vector("list", 5L)
  for (k in 2:6) {
    if (nrow(A) <= k) {
      out[[k-1L]] <- data.frame(k = k, ARI = NA_real_, NMI = NA_real_)
      next
    }
    set.seed(seed + k)
    ar <- tryCatch(kmeans(A, centers = k, nstart = nstart, iter.max = 100)$cluster, error = function(e) NULL)
    set.seed(seed + k)
    br <- tryCatch(kmeans(B, centers = k, nstart = nstart, iter.max = 100)$cluster, error = function(e) NULL)
    out[[k-1L]] <- data.frame(k = k, ARI = if (is.null(ar) || is.null(br)) NA_real_ else ari_score(ar, br), NMI = if (is.null(ar) || is.null(br)) NA_real_ else nmi_score(ar, br))
  }
  do.call(rbind, out)
}

evaluate_set <- function(label, ids, real, pred, group, clinical, de_real_all, de_pred_all, cox_real_all, cox_pred_all, do_cluster = TRUE) {
  R <- real[, ids, drop = FALSE]; P <- pred[, ids, drop = FALSE]
  var_r <- apply(R, 2, var); var_p <- apply(P, 2, var)
  vr <- ifelse(var_r > 0, var_p / var_r, NA_real_)
  mr <- colMeans(R); mp <- colMeans(P)
  rms <- sqrt(colMeans((R-P)^2))
  slope <- vapply(seq_along(ids), function(j) {
    ok <- is.finite(R[, j]) & is.finite(P[, j]); if (sum(ok) < 3 || var(P[ok, j]) == 0) return(NA_real_)
    unname(cov(R[ok, j], P[ok, j]) / var(P[ok, j]))
  }, numeric(1))
  real_de <- de_real_all[match(ids,de_real_all$miRNA),,drop=FALSE]; pred_de <- de_pred_all[match(ids,de_pred_all$miRNA),,drop=FALSE]
  real_c <- cox_real_all[match(ids,cox_real_all$miRNA),,drop=FALSE]; pred_c <- cox_pred_all[match(ids,cox_pred_all$miRNA),,drop=FALSE]
  patient_cor <- vapply(seq_len(nrow(R)), function(i) safe_cor(R[i,],P[i,]), numeric(1))
  corR <- suppressWarnings(cor(R, use = "pairwise.complete.obs")); corP <- suppressWarnings(cor(P, use = "pairwise.complete.obs"))
  upper <- upper.tri(corR)
  cor_abs_diff <- if (any(upper)) mean(abs(corR[upper] - corP[upper]), na.rm = TRUE) else NA_real_
  cor_corr <- if (sum(upper) > 2) safe_cor(corR[upper], corP[upper]) else NA_real_
  de_ok <- is.finite(real_de$effect) & is.finite(pred_de$effect)
  cox_ok <- is.finite(real_c$logHR) & is.finite(pred_c$logHR)
  n_early <- if (is.null(group)) 0L else sum(group == "Early", na.rm = TRUE)
  n_late <- if (is.null(group)) 0L else sum(group == "Late", na.rm = TRUE)
  top10 <- de_overlap(real_de,pred_de,length(ids),10L); top20 <- de_overlap(real_de,pred_de,length(ids),20L); top50 <- de_overlap(real_de,pred_de,length(ids),50L)
  de_slope <- if(sum(de_ok)>2 && var(real_de$effect[de_ok])>0) cov(pred_de$effect[de_ok],real_de$effect[de_ok])/var(real_de$effect[de_ok]) else NA_real_
  summ <- data.frame(subset = label, n_miRNA = length(ids), mean_variance_ratio = mean(vr, na.rm = TRUE), median_variance_ratio = median(vr, na.rm = TRUE), variance_ratio_lt_0_5_fraction = mean(vr < 0.5, na.rm = TRUE), mean_bias = mean(mp-mr, na.rm = TRUE), mean_absolute_bias = mean(abs(mp-mr), na.rm = TRUE), calibration_slope = mean(slope, na.rm = TRUE), mean_RMSE = mean(rms, na.rm = TRUE), DE_effect_correlation = if (sum(de_ok)>2) safe_cor(real_de$effect[de_ok], pred_de$effect[de_ok]) else NA_real_, DE_effect_slope_pred_vs_real = de_slope, DE_direction_agreement = if (sum(de_ok)) mean(sign(real_de$effect[de_ok]) == sign(pred_de$effect[de_ok])) else NA_real_, DE_top10_overlap = top10["overlap"], DE_top10_n_real = top10["real_n"], DE_top10_n_pred = top10["pred_n"], DE_top20_overlap = top20["overlap"], DE_top20_n_real = top20["real_n"], DE_top20_n_pred = top20["pred_n"], DE_top50_overlap = top50["overlap"], DE_top50_n_real = top50["real_n"], DE_top50_n_pred = top50["pred_n"], DE_n_early = n_early, DE_n_late = n_late, DE_BH_lt_0_05_real = sum(p.adjust(real_de$p, "BH") < .05, na.rm=TRUE), DE_BH_lt_0_05_pred = sum(p.adjust(pred_de$p, "BH") < .05, na.rm=TRUE), Cox_logHR_correlation = if (sum(cox_ok)>2) safe_cor(real_c$logHR[cox_ok], pred_c$logHR[cox_ok]) else NA_real_, Cox_logHR_slope = if (sum(cox_ok)>2 && var(pred_c$logHR[cox_ok])>0) cov(real_c$logHR[cox_ok], pred_c$logHR[cox_ok])/var(pred_c$logHR[cox_ok]) else NA_real_, Cox_direction_agreement = if (sum(cox_ok)) mean(sign(real_c$logHR[cox_ok]) == sign(pred_c$logHR[cox_ok])) else NA_real_, Cox_n_valid_real = sum(cox_ok), Cox_events = if (nrow(real_c)) max(real_c$events) else NA_integer_, Cox_censored = if (nrow(real_c)) max(real_c$censored) else NA_integer_, median_patient_correlation = median(patient_cor, na.rm = TRUE), correlation_matrix_correlation = if (sum(upper)>2) cor_corr else NA_real_, correlation_matrix_mean_absolute_difference = cor_abs_diff, check.names = FALSE)
  detail <- data.frame(miRNA = ids, real_mean = mr, predicted_mean = mp, mean_bias = mp-mr, real_variance = var_r, predicted_variance = var_p, variance_ratio = vr, RMSE = rms, calibration_slope = slope, DE_real_trainSD_mean_difference = real_de$effect, DE_pred_trainSD_mean_difference = pred_de$effect, DE_real_p = real_de$p, DE_pred_p = pred_de$p, DE_real_BH = p.adjust(real_de$p, "BH"), DE_pred_BH = p.adjust(pred_de$p, "BH"), Cox_real_logHR_per_training_SD = real_c$logHR, Cox_pred_logHR_per_training_SD = pred_c$logHR, Cox_real_p = real_c$p, Cox_pred_p = pred_c$p, Cox_n = real_c$n, Cox_events = real_c$events, Cox_censored = real_c$censored, stringsAsFactors = FALSE)
  pca <- pca_metrics(R, P)
  cl <- if (do_cluster) cluster_metrics(R, P) else NULL
  list(summary = summ, detail = detail, pca = data.frame(subset=label, t(pca), check.names=FALSE), cluster = if (is.null(cl)) NULL else transform(cl, subset=label))
}

bootstrap_core <- function(real, pred, ids, random_lists, B = 1000L, seed = 261014L) {
  n <- nrow(real)
  if (!length(ids) || n < 5 || B < 1L) return(data.frame())
  core_idx <- match(ids,colnames(real))
  random_idx <- lapply(random_lists,function(s) match(s,colnames(real)))
  if(anyNA(core_idx) || any(vapply(random_idx,anyNA,logical(1)))) stop("Bootstrap feature names do not match matrices")

  # A patient-level correlation across a fixed feature set is unchanged by row
  # resampling; duplicate sampled patients simply duplicate their cached value.
  row_corr <- function(jj) vapply(seq_len(n),function(i) safe_cor(real[i,jj],pred[i,jj]),numeric(1))
  core_row_cor <- row_corr(core_idx)
  random_row_cor <- lapply(random_idx,row_corr)
  patient_median <- function(cor_values, ii) median(cor_values[ii],na.rm=TRUE)

  # Pairwise feature correlations depend on the bootstrap patient multiplicities,
  # but each full 846-feature correlation matrix can be computed once per draw.
  network_abs_diff <- function(ca,cb,jj) {
    if(length(jj)<2L) return(NA_real_)
    delta <- abs(ca[jj,jj,drop=FALSE]-cb[jj,jj,drop=FALSE])
    mean(delta[upper.tri(delta)],na.rm=TRUE)
  }
  network_for_indices <- function(ii) {
    list(real=suppressWarnings(cor(real[ii,,drop=FALSE],use="pairwise.complete.obs")),
         pred=suppressWarnings(cor(pred[ii,,drop=FALSE],use="pairwise.complete.obs")))
  }
  summarize_network <- function(cm) c(core=network_abs_diff(cm$real,cm$pred,core_idx),
    random=vapply(random_idx,function(jj) network_abs_diff(cm$real,cm$pred,jj),numeric(1)))

  # Preserve the exact original per-size seed and the same fixed list of 200 names.
  set.seed(seed); out <- vector("list",B)
  full_rows <- seq_len(n)
  point_core_patient <- patient_median(core_row_cor,full_rows)
  point_random_patient <- vapply(random_row_cor,patient_median,numeric(1),ii=full_rows)
  point_net <- summarize_network(network_for_indices(full_rows))
  point_delta_patient <- point_core_patient-mean(point_random_patient,na.rm=TRUE)
  point_delta_network <- mean(point_net[-1],na.rm=TRUE)-point_net[1]
  for(b in seq_len(B)) {
    ii <- sample.int(n,n,replace=TRUE)
    core_patient <- patient_median(core_row_cor,ii)
    random_patient <- vapply(random_row_cor,patient_median,numeric(1),ii=ii)
    net <- summarize_network(network_for_indices(ii))
    random_patient_mean <- mean(random_patient,na.rm=TRUE)
    random_net_mean <- mean(net[-1],na.rm=TRUE)
    out[[b]] <- data.frame(bootstrap=b,core_median_patient_correlation=core_patient,
      random_mean_median_patient_correlation=random_patient_mean,
      core_correlation_matrix_abs_difference=net[1],
      random_mean_correlation_matrix_abs_difference=random_net_mean,
      Rel_minus_random_median_patient_correlation=core_patient-random_patient_mean,
      random_minus_Rel_correlation_matrix_abs_difference=random_net_mean-net[1],
      point_core_median_patient_correlation=point_core_patient,
      point_random_mean_median_patient_correlation=mean(point_random_patient,na.rm=TRUE),
      point_core_correlation_matrix_abs_difference=point_net[1],
      point_random_mean_correlation_matrix_abs_difference=mean(point_net[-1],na.rm=TRUE),
      point_Rel_minus_random_median_patient_correlation=point_delta_patient,
      point_random_minus_Rel_correlation_matrix_abs_difference=point_delta_network)
  }
  do.call(rbind,out)
}

bootstrap_ci_summary <- function(bt) {
  if(is.null(bt) || !nrow(bt)) return(data.frame())
  boot_metrics <- c("core_median_patient_correlation","core_correlation_matrix_abs_difference","Rel_minus_random_median_patient_correlation","random_minus_Rel_correlation_matrix_abs_difference")
  point_columns <- c("point_core_median_patient_correlation","point_core_correlation_matrix_abs_difference","point_Rel_minus_random_median_patient_correlation","point_random_minus_Rel_correlation_matrix_abs_difference")
  do.call(rbind,lapply(split(bt,bt$requested_size),function(z) do.call(rbind,lapply(boot_metrics,function(nm) {
    x <- z[[nm]]; point_nm <- point_columns[match(nm,boot_metrics)]
    data.frame(requested_size=z$requested_size[1],metric=nm,estimate=z[[point_nm]][1],bootstrap_mean=mean(x,na.rm=TRUE),lower_95=unname(quantile(x,.025,na.rm=TRUE)),upper_95=unname(quantile(x,.975,na.rm=TRUE)),n_boot=sum(is.finite(x)))
  }))))
}

run_e3_bootstrap_only <- function(input_path,output_dir,B=1000L) {
  d <- readRDS(input_path); R <- as_numeric_matrix(d$Y_test_z,"Y_test_z"); P <- as_numeric_matrix(d$pred_test_z,"pred_test_z")
  if(!identical(dimnames(R),dimnames(P))) stop("Y_test_z and pred_test_z dimensions/names must match")
  score <- d$score_table; score <- score[match(colnames(R),score$miRNA),,drop=FALSE]
  if(anyNA(score$miRNA)||any(!is.finite(score$Rel))) stop("Rel scores must cover all test miRNAs")
  rel_order <- score$miRNA[order(-score$Rel,score$miRNA)]
  sizes <- c(33L,74L,177L); rows <- vector("list",length(sizes))
  for(q in seq_along(sizes)) {
    m <- sizes[q]; path <- file.path(output_dir,paste0("E3_random_lists_m",m,".csv"))
    if(!file.exists(path)) stop("Missing completed point-analysis random lists: ",path)
    tab <- read.csv(path,stringsAsFactors=FALSE)
    if(!all(c("random_id","miRNA") %in% names(tab))) stop("Malformed random-list file: ",path)
    tab <- tab[order(tab$random_id),,drop=FALSE]
    by_id <- split(tab$miRNA,tab$random_id)
    if(length(by_id)<200L || any(lengths(by_id[seq_len(200L)])!=m)) stop("Need at least 200 fixed random lists of size ",m)
    random_lists <- lapply(by_id[seq_len(200L)],as.character)
    if(any(vapply(random_lists,function(x) any(!x %in% colnames(R)),logical(1)))) stop("Random-list feature not found in canonical input")
    b <- bootstrap_core(R,P,rel_order[seq_len(m)],random_lists,B=B,seed=20261014L+m)
    b$requested_size <- m; rows[[q]] <- b
  }
  bt <- do.call(rbind,rows); ci <- bootstrap_ci_summary(bt)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  write.csv(bt,file.path(output_dir,paste0("E3_core_patient_bootstrap_fast_B",B,".csv")),row.names=FALSE)
  write.csv(ci,file.path(output_dir,paste0("E3_core_patient_bootstrap_CI_fast_B",B,".csv")),row.names=FALSE)
  saveRDS(list(bootstrap=bt,bootstrap_ci=ci,meta=list(input=normalizePath(input_path,mustWork=TRUE),B=B,seed_by_size=setNames(20261014L+sizes,sizes),random_lists_used=200L,mode="bootstrap_only_fast")),file.path(output_dir,paste0("E3_bootstrap_fast_B",B,".rds")))
  invisible(list(bootstrap=bt,bootstrap_ci=ci))
}

run_e3 <- function(input_path, output_dir, n_random = 1000L, n_boot = 1000L) {
  if (!requireNamespace("survival", quietly=TRUE)) stop("The survival package is required")
  d <- readRDS(input_path)
  required <- c("X_train_z","X_cal_z","X_test_z","Y_train_z","Y_cal_z","Y_test_z","pred_cal_z","pred_test_z","clinical_train","clinical_cal","clinical_test","miRNA_names","gene_names","y_center","y_scale","x_center","x_scale","score_table")
  miss <- setdiff(required,names(d)); if(length(miss)) stop("canonical input missing: ",paste(miss,collapse=", "))
  R <- as_numeric_matrix(d$Y_test_z,"Y_test_z"); P <- as_numeric_matrix(d$pred_test_z,"pred_test_z")
  if(!identical(dim(R),dim(P)) || !identical(dimnames(R),dimnames(P))) stop("Y_test_z and pred_test_z dimensions/names must match")
  if(!identical(colnames(R),as.character(d$miRNA_names))) stop("miRNA_names do not match Y_test_z columns")
  score <- d$score_table
  if(!all(c("miRNA","R2","Pearson","RMSE","Rel") %in% names(score))) stop("score_table requires miRNA,R2,Pearson,RMSE,Rel")
  score <- score[match(colnames(R),score$miRNA),,drop=FALSE]
  if(anyNA(score$miRNA)) stop("score_table must cover all test miRNAs")
  if(any(!is.finite(score$Rel) | !is.finite(score$Pearson) | !is.finite(score$RMSE))) stop("Rel, Pearson, and RMSE must be finite for all test miRNAs")
  clinical <- align_clinical(d$clinical_test,rownames(R),"clinical_test")
  group <- stage_group(clinical,rownames(R)); if(is.null(group)) group <- event_group(clinical,rownames(R))
  group_source <- if(is.null(group)) "not_computable" else attr(group,"source")
  if(is.null(group)) group <- rep(NA_character_,nrow(R))
  de_real_all <- feature_stats(R,group); de_pred_all <- feature_stats(P,group)
  # Cox is fitted once per feature and cached for all subset sizes and random lists.
  cox_real_all <- cox_stats(R,clinical,rownames(R))
  cox_pred_all <- cox_stats(P,clinical,rownames(P))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  set.seed(20261014L)
  choices <- list(Rel=score$miRNA[order(-score$Rel,score$miRNA)], Pearson=score$miRNA[order(-score$Pearson,score$miRNA)], RMSE=score$miRNA[order(score$RMSE,score$miRNA)])
  sizes <- c(33L,74L,177L); all_ids <- colnames(R)
  summary_rows <- pca_rows <- cluster_rows <- rand_rows <- list(); rand_names <- list()
  for(m in sizes) {
    if(m>length(all_ids)) stop("Requested set size exceeds test miRNAs: ",m)
    ids_by <- lapply(choices,function(x) x[seq_len(m)])
    # All and threshold sets are evaluated once; fixed-size comparisons are evaluated at each m.
    if(m == min(sizes)) {
      ids_by$All <- all_ids
      for(th in c(.5,.6,.7)) {
        si <- score$miRNA[is.finite(score$Rel)&score$Rel>=th]
        if(length(si)>0 && length(si)<=length(all_ids)) ids_by[[paste0("Rel_ge_",th)]] <- si
      }
    }
    rnd <- replicate(n_random,sample(all_ids,m,replace=FALSE),simplify=FALSE)
    rand_names[[as.character(m)]] <- rnd[seq_len(min(200L,length(rnd)))]
    write.csv(data.frame(random_id=rep(seq_along(rnd),lengths(rnd)),miRNA=unlist(rnd,use.names=FALSE)),file.path(output_dir,paste0("E3_random_lists_m",m,".csv")),row.names=FALSE)
    names(rnd) <- sprintf("Random_%04d",seq_along(rnd))
    # PCA/clustering are recomputed for all fixed lists; all 1000 random subsets receive the same core PCA/k metrics.
    todo <- c(ids_by,rnd)
    for(nm in names(todo)) {
      z <- evaluate_set(nm,todo[[nm]],R,P,group,clinical,de_real_all,de_pred_all,cox_real_all,cox_pred_all,do_cluster=TRUE)
      metric_size <- length(todo[[nm]])
      z$summary$requested_size <- metric_size; z$summary$group_source <- group_source
      summary_rows[[length(summary_rows)+1L]] <- z$summary
      z$pca$requested_size <- metric_size; pca_rows[[length(pca_rows)+1L]] <- z$pca
      z$cluster$requested_size <- metric_size; cluster_rows[[length(cluster_rows)+1L]] <- z$cluster
      if(grepl("^Random_",nm)) {
        rr <- z$summary; rr$random_id <- as.integer(sub("Random_","",nm)); rand_rows[[length(rand_rows)+1L]] <- rr
      } else {
        write.csv(z$detail,file.path(output_dir,paste0("E3_detail_m",metric_size,"_",gsub("[^A-Za-z0-9_.-]","_",nm),".csv")),row.names=FALSE)
      }
    }
  }
  sm <- do.call(rbind,summary_rows); pc <- do.call(rbind,pca_rows); cl <- do.call(rbind,cluster_rows); rr <- do.call(rbind,rand_rows)
  write.csv(sm,file.path(output_dir,"E3_subset_summary.csv"),row.names=FALSE)
  write.csv(pc,file.path(output_dir,"E3_PCA_metrics.csv"),row.names=FALSE)
  write.csv(cl,file.path(output_dir,"E3_cluster_ARI_NMI.csv"),row.names=FALSE)
  write.csv(rr,file.path(output_dir,"E3_random_subset_distribution.csv"),row.names=FALSE)
  # Core fixed Rel sets: exact sizes 33/74/177 and 200 fixed random lists per size.
  bt <- NULL; boot_ci <- NULL
  bootstrap_status <- if(n_boot>0L) "complete" else "skipped_n_boot_0"
  if(n_boot>0L) {
    boot_rows <- list()
    for(m in sizes) {
      core <- choices$Rel[seq_len(m)]
      b <- bootstrap_core(R,P,core,rand_names[[as.character(m)]],B=n_boot,seed=20261014L+m)
      b$requested_size <- m; boot_rows[[length(boot_rows)+1L]] <- b
    }
    bt <- do.call(rbind,boot_rows); boot_ci <- bootstrap_ci_summary(bt)
    write.csv(bt,file.path(output_dir,paste0("E3_core_patient_bootstrap_B",n_boot,".csv")),row.names=FALSE)
    write.csv(boot_ci,file.path(output_dir,paste0("E3_core_patient_bootstrap_CI_B",n_boot,".csv")),row.names=FALSE)
  } else {
    writeLines("Patient bootstrap was explicitly skipped because n_boot=0. Point analyses and all registered random subsets were still computed; no bootstrap CI is available.",file.path(output_dir,"E3_BOOTSTRAP_SKIPPED.txt"))
  }
  meta <- list(input=normalizePath(input_path,mustWork=TRUE),n_random=n_random,n_boot=n_boot,bootstrap_status=bootstrap_status,seed=20261014L,bootstrap_seed=20261014L,stage_mapping="I/II incl A-C -> Early; III/IV incl A-C -> Late; unknown/stage 0 excluded",group_source=group_source,cox_time_unit="days; all 366 primary times exactly matched to GDC",de_scale="standardized_mean_difference; input miRNA log2(RPM+1) supported by official GDC case comparisons",procrustes="centered normalized patient scores, orthogonal rotation; disparity is mean squared residual",cluster_nstart=50L,created=Sys.time())
  saveRDS(list(summary=sm,pca=pc,cluster=cl,random=rr,bootstrap=bt,bootstrap_ci=boot_ci,meta=meta),file.path(output_dir,"E3_results.rds"))
  writeLines(c("E3 QC",paste0("test patients: ",nrow(R)),paste0("miRNAs: ",ncol(R)),paste0("group source: ",group_source),paste0("group counts: ",paste(capture.output(print(table(group,useNA='ifany'))),collapse=' ')),paste0("sizes: ",paste(sizes,collapse=", ")),paste0("random subsets/size: ",n_random),paste0("bootstrap status: ",bootstrap_status),paste0("cox time unit: days; verified against GDC for all 366 primary cases")),file.path(output_dir,"E3_QC.txt"))
  invisible(list(summary=sm,pca=pc,cluster=cl,random=rr,bootstrap=bt,meta=meta))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly=TRUE)
  n_boot <- as.integer(Sys.getenv("E3_N_BOOT","1000"))
  if(length(args)>=1L && args[1]=="--bootstrap-only") {
    input <- if(length(args)>=2L) args[2] else "new/results/E0/canonical_input.rds"
    out <- if(length(args)>=3L) args[3] else "new/results/E3"
    run_e3_bootstrap_only(input,out,B=n_boot)
  } else {
    input <- if(length(args)>=1L) args[1] else "new/results/E0/canonical_input.rds"
    out <- if(length(args)>=2L) args[2] else "new/results/E3"
    n_random <- as.integer(Sys.getenv("E3_N_RANDOM","1000"))
    run_e3(input,out,n_random=n_random,n_boot=n_boot)
  }
}
