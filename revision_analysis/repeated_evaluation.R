#!/usr/bin/env Rscript
# Post-evaluation for completed E1 repeated-training checkpoints.
# This script never trains models and never reads an in-progress checkpoint.
# Modes: one <seed>, all, watch. All scores are split-local; no cross-seed OOF pooling.

args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)) args[1] else "all"
root <- normalizePath(getwd(),winslash="/",mustWork=TRUE)
e1_dir <- normalizePath(Sys.getenv("E1_INPUT_DIR",file.path(root,"new","results","E1")),winslash="/",mustWork=FALSE)
out_dir <- normalizePath(Sys.getenv("E1_EVALUATION_OUTPUT_DIR",file.path(root,"new","results","E1_evaluation")),winslash="/",mustWork=FALSE)
e0_canonical_path <- Sys.getenv("E1_EVALUATION_E0_CANONICAL",file.path(root,"new","results","E0","canonical_input.rds"))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
seeds_target <- 20261014:20261033
expected_protocol_version <- "trainCV_lambda_default_v2"
evaluation_schema_version <- "e1_eval_schema_v2"
random_subsets_per_size <- as.integer(Sys.getenv("E1_EVALUATION_RANDOM_SUBSETS","200"))
if(is.na(random_subsets_per_size)||random_subsets_per_size<1L) stop("E1_EVALUATION_RANDOM_SUBSETS must be a positive integer")
e0_obj<-readRDS(e0_canonical_path)
e0_names<-e0_obj$miRNA_names
if(is.null(e0_names)||!length(e0_names)) e0_names<-colnames(e0_obj$Y_test_z)
e0_mirna_universe<-unique(as.character(e0_names))
if(!length(e0_mirna_universe)) stop("E0 canonical input does not provide miRNA_names or Y_test_z")
rates <- c(0.039,0.087,0.209)
rate_names <- c("coverage_0039","coverage_0087","coverage_0209")
model_names <- c("ElasticNet","PCA_ElasticNet","PCA_Ridge")
legacy_path <- Sys.getenv("E1_EVALUATION_LEGACY_LIST",file.path(root,"new","results","E0","miRNA_Rel_ge_0.7.csv"))
downstream_path <- file.path(root,"new","code","downstream_analysis.R")
source(downstream_path,local=TRUE)

atomic_save_rds <- function(x,path) {
  tmp <- paste0(path,".tmp")
  saveRDS(x,tmp,compress=TRUE)
  if(file.exists(path)) {
    prior <- paste0(path,".old_",format(Sys.time(),"%Y%m%d%H%M%S"))
    file.rename(path,prior)
  }
  if(!file.rename(tmp,path)) stop("Could not atomically rename checkpoint: ",path)
}
atomic_write_csv <- function(x,path) {
  tmp <- paste0(path,".tmp")
  write.csv(x,tmp,row.names=FALSE,na="")
  if(file.exists(path)) file.remove(path)
  if(!file.rename(tmp,path)) stop("Could not atomically rename CSV: ",path)
}
finite_cor <- function(x,y,method="pearson") {
  ok <- is.finite(x)&is.finite(y)
  if(sum(ok)<3L || sd(x[ok])==0 || sd(y[ok])==0) return(NA_real_)
  suppressWarnings(cor(x[ok],y[ok],method=method))
}
per_feature_metrics <- function(real,pred,scale) {
  real<-as.matrix(real);pred<-as.matrix(pred);scale<-as.numeric(scale)
  if(!identical(dim(real),dim(pred))) stop("real/pred dimensions differ")
  n<-ncol(real);out<-vector("list",n)
  for(j in seq_len(n)) {
    ok<-is.finite(real[,j])&is.finite(pred[,j]);nr<-sum(ok)
    r<-if(nr>=3L) finite_cor(real[ok,j],pred[ok,j]) else NA_real_
    denom<-if(nr>=2L) sum((real[ok,j]-mean(real[ok,j]))^2) else NA_real_
    sse<-if(nr) sum((real[ok,j]-pred[ok,j])^2) else NA_real_
    r2<-if(is.finite(denom)&&denom>0) 1-sse/denom else NA_real_
    rmse<-if(nr) sqrt(mean((real[ok,j]-pred[ok,j])^2)) else NA_real_
    nrmse<-if(nr&&is.finite(scale[j])&&scale[j]>0) rmse/scale[j] else NA_real_
    out[[j]]<-data.frame(miRNA=colnames(real)[j],n_test=nr,test_Pearson=r,test_R2=r2,
      test_RMSE_native=rmse,test_RMSE_trainingSD=nrmse,stringsAsFactors=FALSE)
  }
  do.call(rbind,out)
}
summarize_quality <- function(ids,metrics,seed,model,method,rate,target_n,eligible_n) {
  q<-metrics[match(ids,metrics$miRNA),,drop=FALSE]
  r<-q$test_R2;p<-q$test_Pearson;rn<-q$test_RMSE_native;rz<-q$test_RMSE_trainingSD
  data.frame(seed=seed,model=model,method=method,coverage_rate=rate,n_selected=length(ids),
    target_n=target_n,eligible_n=eligible_n,n_test_r2=sum(is.finite(r)),n_test_pearson=sum(is.finite(p)),
    n_test_rmse=sum(is.finite(rn)),mean_test_R2=if(any(is.finite(r)))mean(r,na.rm=TRUE) else NA_real_,
    median_test_Pearson=if(any(is.finite(p)))median(p,na.rm=TRUE) else NA_real_,
    mean_test_RMSE_native=if(any(is.finite(rn)))mean(rn,na.rm=TRUE) else NA_real_,
    mean_test_RMSE_trainingSD=if(any(is.finite(rz)))mean(rz,na.rm=TRUE) else NA_real_,
    negative_R2_fraction=if(any(is.finite(r)))mean(r<0,na.rm=TRUE) else NA_real_,stringsAsFactors=FALSE)
}
order_ids <- function(tab,score,decreasing=TRUE) {
  ok<-is.finite(score);z<-tab[ok,,drop=FALSE];v<-score[ok]
  if(!nrow(z)) return(character())
  ord<-if(decreasing) order(-v,z$miRNA,method="radix") else order(v,z$miRNA,method="radix")
  z$miRNA[ord]
}
check_complete_seed <- function(seed) {
  path<-file.path(e1_dir,sprintf("seed_%d.rds",seed))
  if(!file.exists(path)) return(NULL)
  x<-tryCatch(readRDS(path),error=function(e)NULL)
  if(is.null(x)||!identical(x$status,"complete")||!identical(as.integer(x$seed),as.integer(seed))) return(NULL)
  if(!identical(x$protocol_version,expected_protocol_version)) return(NULL)
  can<-x$canonical_paths[model_names]
  if(length(can)!=length(model_names)||anyNA(can)||!all(file.exists(can))) return(NULL)
  for(model in model_names) {
    ci<-tryCatch(readRDS(can[[model]]),error=function(e)NULL)
    if(is.null(ci)||!identical(ci$protocol_version,x$protocol_version)||
       !identical(ci$selected_model,model)||!identical(as.integer(ci$split_seed),as.integer(seed))) return(NULL)
  }
  x
}
read_legacy_ids <- function() {
  if(!file.exists(legacy_path)) return(character())
  x<-tryCatch(read.csv(legacy_path,stringsAsFactors=FALSE),error=function(e)NULL)
  if(is.null(x)||!ncol(x)) return(character())
  col<-intersect(c("miRNA","miRNA_name","feature"),names(x))
  if(length(col)) unique(as.character(x[[col[1]]])) else unique(as.character(x[[1]]))
}

evaluate_seed <- function(seed) {
  final_path<-file.path(out_dir,sprintf("seed_%d_evaluation.rds",seed))
  result<-check_complete_seed(seed)
  if(is.null(result)) return(NULL)
  if(file.exists(final_path)) {
    old<-tryCatch(readRDS(final_path),error=function(e)NULL)
    if(!is.null(old)&&identical(old$status,"complete")&&identical(as.integer(old$seed),as.integer(seed))&&
       identical(old$source_protocol_version,result$protocol_version)&&identical(old$evaluation_schema_version,evaluation_schema_version)) return(old)
  }
  message("Evaluating completed E1 seed ",seed)
  t0<-proc.time()[["elapsed"]]
  dir.create(file.path(out_dir,"by_seed"),showWarnings=FALSE,recursive=TRUE)
  seed_dir<-file.path(out_dir,"by_seed",as.character(seed));dir.create(seed_dir,showWarnings=FALSE,recursive=TRUE)
  legacy<-read_legacy_ids()
  per_model<-list();all_assoc<-list();all_scores<-list();all_quality<-list();all_membership<-list();all_threshold<-list();legacy_rows<-list();threshold_membership_rows<-list();performance<-list();hyper<-list()
  canonical<-list()
  for(model in model_names) {
    can_path<-result$canonical_paths[[model]]
    can<-tryCatch(readRDS(can_path),error=function(e)stop("Cannot read canonical checkpoint ",can_path))
    if(!identical(can$selected_model,model)||!identical(as.integer(can$split_seed),as.integer(seed))) stop("Canonical model/seed mismatch: ",can_path)
    if(!all(c("Y_test_native","pred_test_native","Y_test_z","pred_test_z","Y_train_native","score_table","selected_oof")%in%names(can))) stop("Incomplete canonical file: ",can_path)
    if(!identical(colnames(can$Y_test_z),colnames(can$pred_test_z))) stop("Test matrix feature order mismatch for ",model)
    ids<-colnames(can$Y_test_z)
    score<-can$score_table[match(ids,can$score_table$miRNA),,drop=FALSE]
    if(anyNA(score$miRNA)) stop("Calibration score table misses features for ",model)
    oof<-can$selected_oof$normalized_rmse_by_miRNA
    if(!is.data.frame(oof)||!all(c("miRNA","OOF_RMSE_normalized")%in%names(oof))) stop("Missing same-split selected-model OOF baseline for ",model)
    oof<-oof[match(ids,oof$miRNA),,drop=FALSE]
    if(nrow(oof)!=length(ids)) stop("OOF baseline feature length mismatch for ",model)
    train_native<-as.matrix(can$Y_train_native)
    train_mean<-colMeans(train_native,na.rm=TRUE);train_var<-apply(train_native,2,var,na.rm=TRUE)
    test_metrics<-per_feature_metrics(can$Y_test_native,can$pred_test_native,can$y_scale[ids])
    test_metrics_z<-per_feature_metrics(can$Y_test_z,can$pred_test_z,rep(1,length(ids)))
    test_metrics<-test_metrics[match(ids,test_metrics$miRNA),,drop=FALSE]
    # Test-z and native-scale R2/Pearson are mathematically the same, but native RMSE is the primary reported error.
    test_metrics$test_Pearson_z<-test_metrics_z$test_Pearson
    test_metrics$test_R2_z<-test_metrics_z$test_R2
    score$OOF_RMSE_normalized<-oof$OOF_RMSE_normalized
    score$train_native_mean<-train_mean[ids];score$train_native_variance<-train_var[ids]
    score<-merge(score,test_metrics,by="miRNA",sort=FALSE);score<-score[match(ids,score$miRNA),,drop=FALSE]
    score$R2_clipped<-pmax(score$R2,0);score$Pearson_clipped<-pmax(score$Pearson,0)
    score$Rel_recomputed<-0.5*score$R2_clipped+0.5*score$Pearson_clipped
    score$weight_R2_025<-0.25*score$R2_clipped+0.75*score$Pearson_clipped
    score$weight_R2_075<-0.75*score$R2_clipped+0.25*score$Pearson_clipped
    if(any(is.finite(score$Rel)&is.finite(score$Rel_recomputed)&abs(score$Rel-score$Rel_recomputed)>1e-10)) warning("Rel formula mismatch for ",model," seed ",seed)
    eligible_cols<-c("Rel","R2","Pearson","RMSE","OOF_RMSE_normalized","train_native_mean","train_native_variance")
    common_eval<-Reduce(`&`,lapply(eligible_cols,function(nm)is.finite(score[[nm]])))
    common_ids<-score$miRNA[common_eval];n_eval<-length(common_ids)
    sizes<-as.integer(floor(n_eval*rates+0.5));names(sizes)<-rate_names
    if(!n_eval||any(sizes<1L)) stop("No common eligible features at registered coverages for seed/model ",seed,"/",model)
    methods<-list(
      Rel=list(score=score$Rel,direction="high"),
      raw_R2=list(score=score$R2,direction="high"),clipped_R2=list(score=score$R2_clipped,direction="high"),
      raw_Pearson=list(score=score$Pearson,direction="high"),clipped_Pearson=list(score=score$Pearson_clipped,direction="high"),
      weight_R2_025_Pearson_075=list(score=score$weight_R2_025,direction="high"),
      weight_R2_075_Pearson_025=list(score=score$weight_R2_075,direction="high"),
      calibration_RMSE=list(score=score$RMSE,direction="low"),
      same_split_training_OOF_RMSE=list(score=score$OOF_RMSE_normalized,direction="low"),
      training_native_mean=list(score=score$train_native_mean,direction="high"),
      training_native_variance=list(score=score$train_native_variance,direction="high"))
    list_for_model<-list();membership<-list();quality<-list();assoc<-list()
    test_error<-score$test_RMSE_trainingSD
    test_error_native<-score$test_RMSE_native
    for(method in names(methods)) {
      sc<-methods[[method]]$score;dec<-methods[[method]]$direction=="high"
      score_quality<-if(dec) sc else -sc
      ord<-order_ids(score,sc,decreasing=dec);ord<-ord[ord%in%common_ids]
      # score list order is stable: score then lexical miRNA tie break.
      for(q in seq_along(rates)) {
        k<-sizes[q];chosen<-head(ord,k);rate_name<-rate_names[q]
        list_for_model[[paste(method,rate_name,sep="__")]]<-chosen
        quality[[length(quality)+1L]]<-summarize_quality(chosen,score,seed,model,method,rates[q],k,n_eval)
        for(ii in seq_along(score$miRNA)) {
          mir<-score$miRNA[ii];is_eval<-common_eval[ii]
          membership[[length(membership)+1L]]<-data.frame(seed=seed,model=model,method=method,coverage_rate=rates[q],target_n=k,eligible_n=n_eval,miRNA=mir,score=sc[ii],evaluable=is_eval,selected=is_eval&&mir%in%chosen,stringsAsFactors=FALSE)
        }
      }
      assoc[[length(assoc)+1L]]<-data.frame(seed=seed,model=model,method=method,score_direction=methods[[method]]$direction,
        n_score_test_RMSE_trainingSD=sum(is.finite(sc)&is.finite(test_error)),
        spearman_raw_score_vs_test_RMSE_trainingSD=finite_cor(sc,test_error,"spearman"),
        spearman_quality_score_vs_test_RMSE_trainingSD=finite_cor(score_quality,test_error,"spearman"),
        n_score_test_RMSE_native=sum(is.finite(sc)&is.finite(test_error_native)),
        spearman_raw_score_vs_test_RMSE_native=finite_cor(sc,test_error_native,"spearman"),
        spearman_quality_score_vs_test_RMSE_native=finite_cor(score_quality,test_error_native,"spearman"),
        n_score_test_R2=sum(is.finite(sc)&is.finite(score$test_R2)),
        spearman_raw_score_vs_test_R2=finite_cor(sc,score$test_R2,"spearman"),
        spearman_quality_score_vs_test_R2=finite_cor(score_quality,score$test_R2,"spearman"),
        n_score_test_Pearson=sum(is.finite(sc)&is.finite(score$test_Pearson)),
        spearman_raw_score_vs_test_Pearson=finite_cor(sc,score$test_Pearson,"spearman"),
        spearman_quality_score_vs_test_Pearson=finite_cor(score_quality,score$test_Pearson,"spearman"),stringsAsFactors=FALSE)
    }
    # Actual Rel thresholds remain separate from fixed-coverage lists.
    for(th in c(.5,.6,.7)) {
      chosen<-score$miRNA[is.finite(score$Rel)&score$Rel>=th]
      n_rel_eval<-sum(is.finite(score$Rel))
      q<-summarize_quality(chosen,score,seed,model,paste0("Rel_ge_",th),NA_real_,length(chosen),n_rel_eval)
      quality[[length(quality)+1L]]<-q
      all_threshold[[length(all_threshold)+1L]]<-data.frame(seed=seed,model=model,selection_kind="Rel_threshold",threshold=th,n_selected=length(chosen),n_eligible=n_rel_eval,actual_coverage=if(n_rel_eval)length(chosen)/n_rel_eval else NA_real_,stringsAsFactors=FALSE)
      list_for_model[[paste0("Rel_ge_",th)]]<-chosen
      for(mir in e0_mirna_universe) {
        idx<-match(mir,score$miRNA);rel_eval<-!is.na(idx)&&is.finite(score$Rel[idx]);rel_value<-if(rel_eval)score$Rel[idx]else NA_real_
        threshold_membership_rows[[length(threshold_membership_rows)+1L]]<-data.frame(seed=seed,model=model,selection_kind="Rel_threshold",threshold=th,miRNA=mir,score=rel_value,evaluable=rel_eval,selected=rel_eval&&rel_value>=th,stringsAsFactors=FALSE)
      }
    }
    # Operational labels are defined from test measurements; they are used only to calculate recall/FNR.
    for(label_r in c(.3,.5,.7)) {
      positive<-is.finite(score$test_R2)&score$test_R2>0&is.finite(score$test_Pearson)&score$test_Pearson>=label_r
      n_positive<-sum(positive)
      for(method in names(methods)) for(q in seq_along(rates)) {
        chosen<-list_for_model[[paste(method,rate_names[q],sep="__")]]
        hit<-sum(chosen%in%score$miRNA[positive]);recall<-if(n_positive)hit/n_positive else NA_real_
        all_scores[[length(all_scores)+1L]]<-data.frame(seed=seed,model=model,method=method,coverage_rate=rates[q],target_n=sizes[q],
          operating_r=label_r,n_positive=n_positive,n_selected=length(chosen),n_positive_selected=hit,recall=recall,FNR=if(n_positive)1-recall else NA_real_,stringsAsFactors=FALSE)
      }
      for(th in c(.5,.6,.7)) {
        chosen<-list_for_model[[paste0("Rel_ge_",th)]]
        hit<-sum(chosen%in%score$miRNA[positive]);recall<-if(n_positive)hit/n_positive else NA_real_
        all_scores[[length(all_scores)+1L]]<-data.frame(seed=seed,model=model,method=paste0("Rel_ge_",th),coverage_rate=length(chosen)/n_eval,target_n=length(chosen),
          operating_r=label_r,n_positive=n_positive,n_selected=length(chosen),n_positive_selected=hit,recall=recall,FNR=if(n_positive)1-recall else NA_real_,stringsAsFactors=FALSE)
      }
    }
    # Preserve original 33-item membership and report evaluation/selection per item and method.
    if(length(legacy)) for(mir in legacy) {
      idx<-match(mir,score$miRNA);is_eval<-!is.na(idx)&&common_eval[idx]
      rel_eval<-!is.na(idx)&&is.finite(score$Rel[idx])
      row<-data.frame(seed=seed,model=model,miRNA=mir,evaluable=is_eval,Rel_evaluable=rel_eval,
        test_R2=if(is.na(idx))NA_real_ else score$test_R2[idx],test_Pearson=if(is.na(idx))NA_real_ else score$test_Pearson[idx],
        test_RMSE_native=if(is.na(idx))NA_real_ else score$test_RMSE_native[idx],Rel=if(is.na(idx))NA_real_ else score$Rel[idx],
        Rel_ge_0_5=rel_eval&&score$Rel[idx]>=.5,Rel_ge_0_6=rel_eval&&score$Rel[idx]>=.6,Rel_ge_0_7=rel_eval&&score$Rel[idx]>=.7,
        Rel_top_0039=mir%in%list_for_model[[paste("Rel",rate_names[1],sep="__")]],
        rawR2_top_0039=mir%in%list_for_model[[paste("raw_R2",rate_names[1],sep="__")]],
        rawPearson_top_0039=mir%in%list_for_model[[paste("raw_Pearson",rate_names[1],sep="__")]],
        calibrationRMSE_top_0039=mir%in%list_for_model[[paste("calibration_RMSE",rate_names[1],sep="__")]],
        same_split_OOF_top_0039=mir%in%list_for_model[[paste("same_split_training_OOF_RMSE",rate_names[1],sep="__")]],stringsAsFactors=FALSE)
      legacy_rows[[length(legacy_rows)+1L]]<-row
    }
    selected_row<-result$selected[result$selected$model==model,,drop=FALSE]
    selected_row$seed<-seed;selected_row$winner_flag<-selected_row$model==result$selected_model;hyper[[length(hyper)+1L]]<-selected_row

    per_model[[model]]<-list(score=score,eligible_ids=common_ids,sizes=sizes,lists=list_for_model,
      quality=do.call(rbind,quality),associations=do.call(rbind,assoc),
      legacy=if(length(legacy_rows)) do.call(rbind,legacy_rows[vapply(legacy_rows,function(x)x$model[1]==model,logical(1))]) else data.frame(),
      membership=if(length(membership)) do.call(rbind,membership) else data.frame())
    all_assoc[[model]]<-do.call(rbind,assoc);all_quality[[model]]<-do.call(rbind,quality);all_membership[[model]]<-do.call(rbind,membership)
    canonical[[model]]<-can
  }
  # Assemble E2 tables across model-specific canonical files.
  quality_tab<-do.call(rbind,lapply(per_model,`[[`,"quality"))
  assoc_tab<-do.call(rbind,lapply(per_model,`[[`,"associations"))
  membership_tab<-do.call(rbind,lapply(per_model,`[[`,"membership"))
  # Recreate FNR rows once in a direct, auditable pass from model score tables.
  fnr_rows<-list()
  for(model in model_names) {
    z<-per_model[[model]]$score;lists<-per_model[[model]]$lists
    for(label_r in c(.3,.5,.7)) {
      pos<-is.finite(z$test_R2)&z$test_R2>0&is.finite(z$test_Pearson)&z$test_Pearson>=label_r;npos<-sum(pos)
      for(method in names(per_model[[model]]$lists)[grepl("__coverage_",names(per_model[[model]]$lists))]) {
        ss<-strsplit(method,"__",fixed=TRUE)[[1]];chosen<-lists[[method]];hit<-sum(chosen%in%z$miRNA[pos]);rec<-if(npos)hit/npos else NA_real_
        fnr_rows[[length(fnr_rows)+1L]]<-data.frame(seed=seed,model=model,method=ss[1],coverage=ss[2],operating_r=label_r,n_positive=npos,n_selected=length(chosen),n_positive_selected=hit,recall=rec,FNR=if(npos)1-rec else NA_real_)
      }
      for(th in c(.5,.6,.7)) {
        method<-paste0("Rel_ge_",th);chosen<-lists[[method]];hit<-sum(chosen%in%z$miRNA[pos]);rec<-if(npos)hit/npos else NA_real_
        fnr_rows[[length(fnr_rows)+1L]]<-data.frame(seed=seed,model=model,method=method,coverage=length(chosen)/length(per_model[[model]]$eligible_ids),operating_r=label_r,n_positive=npos,n_selected=length(chosen),n_positive_selected=hit,recall=rec,FNR=if(npos)1-rec else NA_real_)
      }
    }
  }
  all_scores_tab<-do.call(rbind,fnr_rows)
  legacy_tab<-if(length(legacy_rows))do.call(rbind,legacy_rows)else data.frame()
  performance_tab<-result$test_metrics[result$test_metrics$model%in%c(model_names,"Mean"),,drop=FALSE]
  performance_tab$seed<-seed;performance_tab$training_selected_model<-result$selected_model;performance_tab$training_selected_model_flag<-performance_tab$model==result$selected_model
  performance_tab$rank_median_pearson<-rank(-performance_tab$median_pearson,ties.method="min",na.last="keep")
  performance_tab$rank_mean_r2<-rank(-performance_tab$mean_r2,ties.method="min",na.last="keep")
  performance_tab$rank_normalized_rmse<-rank(performance_tab$normalized_rmse,ties.method="min",na.last="keep")
  performance_tab$selected_config<-NA_integer_;performance_tab$selected_k<-NA_integer_;performance_tab$selected_alpha<-NA_real_;performance_tab$selected_lambda_path_index<-NA_integer_;performance_tab$selected_lambda_ratio<-NA_real_
  for(model in model_names) {h<-result$selected[result$selected$model==model,,drop=FALSE];ii<-match(model,performance_tab$model);if(nrow(h)){performance_tab[ii,c("selected_config","selected_k","selected_alpha","selected_lambda_path_index","selected_lambda_ratio")]<-h[1,c("config","k","alpha","lambda_path_index","lambda_ratio")]}}
  hyper_tab<-do.call(rbind,hyper)
  # E3 repeated-split analysis is limited to the training-OOF-selected winner.
  winner<-result$selected_model
  can<-canonical[[winner]];R<-as.matrix(can$Y_test_z);P<-as.matrix(can$pred_test_z)
  clinical<-align_clinical(can$clinical_test,rownames(R),"clinical_test")
  grp<-stage_group(clinical,rownames(R));if(is.null(grp))grp<-event_group(clinical,rownames(R));if(is.null(grp))grp<-rep(NA_character_,nrow(R))
  wid<-per_model[[winner]]$score;univ<-per_model[[winner]]$eligible_ids;n_eval<-length(univ);sizes<-per_model[[winner]]$sizes
  de_r<-feature_stats(R,grp);de_p<-feature_stats(P,grp);ys<-can$y_scale;ys<-ys[colnames(R)]
  cxr<-cox_stats(R,clinical,rownames(R));cxp<-cox_stats(P,clinical,rownames(P))
  e3_selectors<-list(Rel=wid$Rel,Pearson=wid$Pearson,calRMSE=wid$RMSE,training_OOF_RMSE=wid$OOF_RMSE_normalized)
  e3_out<-list(summary=list(),pca=list(),cluster=list(),random_lists=list(),threshold=list())
  set.seed(as.integer(seed)+700000L)
  for(q in seq_along(rates)) {
    k<-sizes[q];lists<-list()
    for(nm in names(e3_selectors)) {
      sc<-e3_selectors[[nm]];ord<-order_ids(wid,sc,decreasing=nm%in%c("Rel","Pearson"));lists[[nm]]<-head(ord[ord%in%univ],k)
    }
    random<-replicate(random_subsets_per_size,sample(univ,k,replace=FALSE),simplify=FALSE)
    names(random)<-sprintf("Random_%03d",seq_along(random));lists<-c(lists,random)
    atomic_write_csv(data.frame(random_id=rep(seq_along(random),lengths(random)),miRNA=unlist(random,use.names=FALSE)),file.path(seed_dir,paste0("E3_random_lists_m",k,".csv")))
    e3_out$random_lists[[as.character(k)]]<-random
    for(nm in names(lists)) {
      ev<-evaluate_set(nm,lists[[nm]],R,P,grp,clinical,de_r,de_p,cxr,cxp,do_cluster=TRUE)
      ev$summary$seed<-seed;ev$summary$training_selected_model<-winner;ev$summary$requested_coverage<-rates[q]
      ev$pca$seed<-seed;ev$pca$training_selected_model<-winner;ev$pca$requested_coverage<-rates[q]
      ev$cluster$seed<-seed;ev$cluster$training_selected_model<-winner;ev$cluster$requested_coverage<-rates[q]
      ev$summary$random_id<-if(grepl("^Random_",nm))as.integer(sub("Random_","",nm)) else NA_integer_
      e3_out$summary[[length(e3_out$summary)+1L]]<-ev$summary;e3_out$pca[[length(e3_out$pca)+1L]]<-ev$pca;e3_out$cluster[[length(e3_out$cluster)+1L]]<-ev$cluster
    }
  }
  for(th in c(.5,.6,.7)) {
    ids_th<-wid$miRNA[is.finite(wid$Rel)&wid$Rel>=th]
    if(length(ids_th)) {
      ev<-evaluate_set(paste0("Rel_ge_",th),ids_th,R,P,grp,clinical,de_r,de_p,cxr,cxp,do_cluster=TRUE)
      ev$summary$seed<-seed;ev$summary$training_selected_model<-winner;ev$summary$requested_coverage<-NA_real_
      e3_out$threshold[[length(e3_out$threshold)+1L]]<-ev$summary
    }
  }
  e3_tab<-do.call(rbind,e3_out$summary);e3_pca<-do.call(rbind,e3_out$pca);e3_cluster<-do.call(rbind,e3_out$cluster);e3_threshold<-if(length(e3_out$threshold))do.call(rbind,e3_out$threshold)else data.frame()
  threshold_membership_tab<-do.call(rbind,threshold_membership_rows)
  # Save compact CSV artifacts before the completion checkpoint; retries can safely overwrite own files.
  atomic_write_csv(performance_tab,file.path(seed_dir,"model_performance_and_test_ranks.csv"))
  atomic_write_csv(hyper_tab,file.path(seed_dir,"selected_hyperparameters.csv"))
  atomic_write_csv(quality_tab,file.path(seed_dir,"selection_quality.csv"))
  atomic_write_csv(assoc_tab,file.path(seed_dir,"score_test_associations.csv"))
  atomic_write_csv(all_scores_tab,file.path(seed_dir,"operational_FNR.csv"))
  atomic_write_csv(legacy_tab,file.path(seed_dir,"original_33_membership.csv"))
  atomic_write_csv(membership_tab,file.path(seed_dir,"selection_membership_long.csv"))
  atomic_write_csv(threshold_membership_tab,file.path(seed_dir,"Rel_threshold_membership_long.csv"))
  atomic_write_csv(do.call(rbind,all_threshold),file.path(seed_dir,"Rel_threshold_counts.csv"))
  atomic_write_csv(e3_tab,file.path(seed_dir,"E3_winner_equal_size_summary.csv"))
  atomic_write_csv(e3_pca,file.path(seed_dir,"E3_winner_PCA.csv"))
  atomic_write_csv(e3_cluster,file.path(seed_dir,"E3_winner_ARI_NMI.csv"))
  atomic_write_csv(e3_threshold,file.path(seed_dir,"E3_winner_Rel_thresholds.csv"))
  for(model in model_names) atomic_write_csv(per_model[[model]]$score,file.path(seed_dir,paste0("miRNA_score_test_metrics_",model,".csv")))
  compact<-list(status="complete",evaluation_schema_version=evaluation_schema_version,seed=seed,source_protocol_version=result$protocol_version,selected_model=winner,split_counts=result$split_counts,
    performance=performance_tab,hyperparameters=hyper_tab,quality=quality_tab,associations=assoc_tab,FNR=all_scores_tab,
    legacy_membership=legacy_tab,membership=membership_tab,threshold_membership=threshold_membership_tab,threshold_counts=do.call(rbind,all_threshold),
    E3_summary=e3_tab,E3_PCA=e3_pca,E3_cluster=e3_cluster,E3_threshold=e3_threshold,
    E3_random_lists=e3_out$random_lists,eligible_n_by_model=setNames(vapply(per_model,function(x)length(x$eligible_ids),integer(1)),model_names),
    elapsed_seconds=proc.time()[["elapsed"]]-t0,metadata=list(coverage_rates=rates,coverage_sizes_by_model=lapply(per_model,`[[`,"sizes"),
      E3_n_boot=0L,random_subsets_per_size=random_subsets_per_size,OOF_baseline="same split, same selected model, training OOF normalized RMSE by miRNA",tie_break="lexical miRNA name",created=Sys.time()))
  atomic_save_rds(compact,final_path)
  message("Evaluation complete for seed ",seed,"; winner=",winner,"; elapsed ",round(compact$elapsed_seconds/60,1)," min")
  compact
}

aggregate_results <- function() {
  eval_files<-list.files(out_dir,pattern="^seed_[0-9]+_evaluation\\.rds$",full.names=TRUE)
  objs<-lapply(eval_files,function(f)tryCatch(readRDS(f),error=function(e)NULL));objs<-Filter(function(x)!is.null(x)&&identical(x$status,"complete")&&identical(x$source_protocol_version,expected_protocol_version)&&identical(x$evaluation_schema_version,evaluation_schema_version),objs)
  done<-sort(vapply(objs,`[[`,integer(1),"seed"))
  if(!length(objs)) {
    state<-data.frame(protocol_version=expected_protocol_version,evaluation_schema_version=evaluation_schema_version,target_seeds=length(seeds_target),completed_seed_evaluations=0L,completed_seeds="",status="waiting_for_protocol_v2_complete_checkpoints",last_updated=as.character(Sys.time()))
    atomic_write_csv(state,file.path(out_dir,"E1_evaluation_progress.csv"))
    message("E1 evaluation: 0/20 complete; waiting for protocol ",expected_protocol_version," checkpoints")
    return(invisible(integer()))
  }
  perf<-do.call(rbind,lapply(objs,`[[`,"performance"));hyp<-do.call(rbind,lapply(objs,`[[`,"hyperparameters"));qual<-do.call(rbind,lapply(objs,`[[`,"quality"));mem<-do.call(rbind,lapply(objs,`[[`,"membership"));thmem<-do.call(rbind,lapply(objs,`[[`,"threshold_membership"));legacy<-do.call(rbind,lapply(objs,`[[`,"legacy_membership"));e3<-do.call(rbind,lapply(objs,`[[`,"E3_summary"))
  atomic_write_csv(perf,file.path(out_dir,"E1_model_performance_by_seed.csv"))
  atomic_write_csv(hyp,file.path(out_dir,"E1_selected_hyperparameters_by_seed.csv"))
  atomic_write_csv(qual,file.path(out_dir,"E1_selection_quality_by_seed.csv"))
  atomic_write_csv(e3,file.path(out_dir,"E3_winner_by_seed.csv"))
  hyp_fields<-c("model","config","k","alpha","lambda_path_index","lambda_ratio")
  hyp_key<-do.call(paste,c(lapply(hyp[hyp_fields],function(v)ifelse(is.na(v),"<NA>",as.character(v))),sep="\r"))
  hyp_group<-unique(hyp_key);hyp_freq<-hyp[match(hyp_group,hyp_key),hyp_fields,drop=FALSE]
  hyp_freq$n_seeds<-vapply(hyp_group,function(k)sum(hyp_key==k),integer(1))
  hyp_freq$selection_frequency<-hyp_freq$n_seeds/length(done)
  atomic_write_csv(hyp_freq,file.path(out_dir,"E1_hyperparameter_selection_frequency.csv"))
  freq<-aggregate(cbind(evaluated=as.integer(mem$evaluable),selected=as.integer(mem$selected))~model+method+coverage_rate+miRNA,mem,sum)
  colnames(freq)[colnames(freq)=="evaluated"]<-"n_evaluable_splits";colnames(freq)[colnames(freq)=="selected"]<-"n_selected_splits"
  freq$selection_frequency<-freq$n_selected_splits/freq$n_evaluable_splits
  atomic_write_csv(freq,file.path(out_dir,"E1_miRNA_selection_frequency.csv"))
  thfreq<-aggregate(cbind(evaluable=as.integer(thmem$evaluable),selected=as.integer(thmem$selected))~model+threshold+miRNA,thmem,sum)
  colnames(thfreq)[colnames(thfreq)=="evaluable"]<-"n_evaluable_splits";colnames(thfreq)[colnames(thfreq)=="selected"]<-"n_selected_splits"
  thfreq$n_recorded_splits<-length(done);thfreq$n_missing_splits<-length(done)-thfreq$n_evaluable_splits
  thfreq$selection_frequency<-ifelse(thfreq$n_evaluable_splits>0,thfreq$n_selected_splits/thfreq$n_evaluable_splits,NA_real_)
  atomic_write_csv(thfreq,file.path(out_dir,"E1_Rel_threshold_selection_frequency.csv"))
  orig_keys<-unique(legacy[c("model","miRNA")]);orig_rows<-vector("list",nrow(orig_keys))
  orig_select_cols<-c("Rel_top_0039","rawR2_top_0039","rawPearson_top_0039","calibrationRMSE_top_0039","same_split_OOF_top_0039")
  for(i in seq_len(nrow(orig_keys))) {
    z<-legacy[legacy$model==orig_keys$model[i]&legacy$miRNA==orig_keys$miRNA[i],,drop=FALSE]
    ne<-sum(z$evaluable);nr<-sum(z$Rel_evaluable)
    vals<-list(model=orig_keys$model[i],miRNA=orig_keys$miRNA[i],n_recorded_splits=length(done),n_evaluable_splits=ne,n_missing_splits=length(done)-ne,n_Rel_evaluable_splits=nr,n_Rel_missing_splits=length(done)-nr)
    for(nm in orig_select_cols) {nn<-sum(z[[nm]]&z$evaluable);vals[[paste0("n_",nm)]]<-nn;vals[[paste0(nm,"_frequency")]]<-if(ne>0)nn/ne else NA_real_}
    for(th in c(.5,.6,.7)) {nm<-paste0("Rel_ge_0_",as.integer(th*10));nn<-sum(z[[nm]]&z$Rel_evaluable);vals[[paste0("n_",nm)]]<-nn;vals[[paste0(nm,"_frequency")]]<-if(nr>0)nn/nr else NA_real_}
    orig_rows[[i]]<-as.data.frame(vals,stringsAsFactors=FALSE)
  }
  atomic_write_csv(do.call(rbind,orig_rows),file.path(out_dir,"E1_original_33_selection_frequency.csv"))
  # Threshold-list overlap uses the actual Rel cutoffs and keeps threshold as its own key.
  thjacc<-list();thkeys<-unique(thmem[c("model","threshold")])
  for(i in seq_len(nrow(thkeys))) {
    z<-thkeys[i,];ss<-sort(unique(thmem$seed[thmem$model==z$model&thmem$threshold==z$threshold]))
    if(length(ss)<2L) next
    for(a in seq_len(length(ss)-1L)) for(b in (a+1L):length(ss)) {
      x<-thmem$miRNA[thmem$seed==ss[a]&thmem$model==z$model&thmem$threshold==z$threshold&thmem$selected]
      y<-thmem$miRNA[thmem$seed==ss[b]&thmem$model==z$model&thmem$threshold==z$threshold&thmem$selected]
      thjacc[[length(thjacc)+1L]]<-data.frame(model=z$model,selection_kind="Rel_threshold",threshold=z$threshold,seed_a=ss[a],seed_b=ss[b],n_a=length(x),n_b=length(y),intersection=length(intersect(x,y)),union=length(union(x,y)),Jaccard=if(length(union(x,y)))length(intersect(x,y))/length(union(x,y))else NA_real_)
    }
  }
  if(length(thjacc))atomic_write_csv(do.call(rbind,thjacc),file.path(out_dir,"E1_Rel_threshold_pairwise_Jaccard.csv"))
  # Pairwise Jaccard is descriptive across overlapping repeated splits, never treated as independent trials.
  jacc<-list();by_method<-unique(mem[c("model","method","coverage_rate","seed")])
  keys<-unique(by_method[c("model","method","coverage_rate")])
  for(i in seq_len(nrow(keys))) {
    z<-keys[i,];ss<-sort(unique(mem$seed[mem$model==z$model&mem$method==z$method&mem$coverage_rate==z$coverage_rate]))
    if(length(ss)<2L) next
    for(a in seq_len(length(ss)-1L)) for(b in (a+1L):length(ss)) {
      x<-mem$miRNA[mem$seed==ss[a]&mem$model==z$model&mem$method==z$method&mem$coverage_rate==z$coverage_rate&mem$selected]
      y<-mem$miRNA[mem$seed==ss[b]&mem$model==z$model&mem$method==z$method&mem$coverage_rate==z$coverage_rate&mem$selected]
      jacc[[length(jacc)+1L]]<-data.frame(model=z$model,method=z$method,coverage_rate=z$coverage_rate,seed_a=ss[a],seed_b=ss[b],n_a=length(x),n_b=length(y),intersection=length(intersect(x,y)),union=length(union(x,y)),Jaccard=if(length(union(x,y)))length(intersect(x,y))/length(union(x,y)) else NA_real_)
    }
  }
  if(length(jacc)) atomic_write_csv(do.call(rbind,jacc),file.path(out_dir,"E1_pairwise_Jaccard.csv"))
  winner_freq<-as.data.frame(table(vapply(objs,`[[`,character(1),"selected_model")),stringsAsFactors=FALSE);names(winner_freq)<-c("selected_model","n_seeds")
  winner_freq$frequency<-winner_freq$n_seeds/length(objs);atomic_write_csv(winner_freq,file.path(out_dir,"E1_winner_frequency_partial.csv"))
  # Summaries remain explicitly partial until all 20 pre-registered seeds have complete evaluation checkpoints.
  state<-data.frame(protocol_version=expected_protocol_version,evaluation_schema_version=evaluation_schema_version,target_seeds=length(seeds_target),completed_seed_evaluations=length(done),completed_seeds=paste(done,collapse=";"),
    status=if(length(done)==length(seeds_target))"all_20_complete" else "partial_no_final_20_seed_conclusion",last_updated=as.character(Sys.time()))
  atomic_write_csv(state,file.path(out_dir,"E1_evaluation_progress.csv"))
  if(length(done)==length(seeds_target)) {
    atomic_write_csv(aggregate(perf[c("median_pearson","mean_r2","normalized_rmse","overall_rmse")],list(model=perf$model),function(x)c(mean=mean(x,na.rm=TRUE),median=median(x,na.rm=TRUE),q05=quantile(x,.05,na.rm=TRUE),q95=quantile(x,.95,na.rm=TRUE))),file.path(out_dir,"E1_final_model_performance_distribution.csv"))
  }
  invisible(done)
}

process_available <- function(requested=seeds_target) {
  done<-integer()
  for(seed in requested) {
    if(!is.null(evaluate_seed(seed))) done<-c(done,seed)
  }
  aggregate_results();done
}

if(mode=="watch") {
  repeat {
    done<-process_available()
    if(length(done)>=length(seeds_target) || file.exists(file.path(out_dir,"E1_evaluation_progress.csv")) && read.csv(file.path(out_dir,"E1_evaluation_progress.csv"))$status[1]=="all_20_complete") break
    message("Completed-evaluation checkpoints processed: ",length(done)," this scan; waiting 60 seconds for complete E1 checkpoints")
    Sys.sleep(60)
  }
} else if(mode=="all") {
  process_available()
} else {
  requested<-suppressWarnings(as.integer(strsplit(mode,",",fixed=TRUE)[[1]]))
  if(anyNA(requested)||any(!requested%in%seeds_target)) stop("Usage: Rscript new/code/repeated_evaluation.R [all|watch|seed[,seed...]]")
  process_available(requested)
}

