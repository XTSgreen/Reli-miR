#!/usr/bin/env Rscript
# Reproduce the patient-paired bootstrap contrasts in E3_paired_comparison_bootstrap_CI.csv.
# This script reuses the fixed canonical test matrices and the first 200 fixed
# random lists at each size; it does not regenerate model fits or feature lists.

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(getwd(),winslash="/",mustWork=TRUE)
input_path <- if(length(args)>=1L) args[[1L]] else file.path(root,"new","results","E0","canonical_input.rds")
list_dir <- if(length(args)>=2L) args[[2L]] else file.path(root,"new","results","E3")
output_path <- if(length(args)>=3L) args[[3L]] else file.path(list_dir,"E3_paired_comparison_bootstrap_CI.csv")
B <- if(length(args)>=4L) as.integer(args[[4L]]) else 1000L
if(is.na(B)||B<1L) stop("B must be a positive integer")

safe_cor <- function(x,y) {
  ok <- is.finite(x)&is.finite(y)
  if(sum(ok)<3L||sd(x[ok])==0||sd(y[ok])==0) return(NA_real_)
  suppressWarnings(cor(x[ok],y[ok]))
}
patient_cor_by_row <- function(real,pred,idx) {
  vapply(seq_len(nrow(real)),function(i)safe_cor(real[i,idx],pred[i,idx]),numeric(1))
}
patient_median <- function(v,ii) median(v[ii],na.rm=TRUE)
network_pairs <- function(id_sets,feature_names) {
  lapply(id_sets,function(ids) {
    ix<-match(ids,feature_names)
    if(anyNA(ix)||length(ix)<2L) stop("Feature list cannot be aligned or has fewer than two items")
    ij<-which(upper.tri(matrix(FALSE,length(ix),length(ix))),arr.ind=TRUE)
    list(i=ix[ij[,1]],j=ix[ij[,2]])
  })
}
network_abs_diff <- function(cr,cp,pairs) mean(abs(cr[cbind(pairs$i,pairs$j)]-cp[cbind(pairs$i,pairs$j)]),na.rm=TRUE)
network_values <- function(cr,cp,pair_sets) vapply(pair_sets,function(z)network_abs_diff(cr,cp,z),numeric(1))

d<-readRDS(input_path)
real<-as.matrix(d$Y_test_z);pred<-as.matrix(d$pred_test_z)
if(!identical(dimnames(real),dimnames(pred))) stop("Y_test_z and pred_test_z must have identical dimnames")
score<-d$score_table
if(!all(c("miRNA","Rel","Pearson","RMSE")%in%names(score))) stop("score_table needs miRNA, Rel, Pearson, RMSE")
score<-score[match(colnames(real),score$miRNA),,drop=FALSE]
if(anyNA(score$miRNA)||any(!is.finite(score$Rel)|!is.finite(score$Pearson)|!is.finite(score$RMSE))) stop("Required scores must cover all test features")
sets_by_size<-list();random_by_size<-list()
sizes<-c(33L,74L,177L)
reference_path<-file.path(list_dir,"E3_paired_comparison_bootstrap_CI.csv")
reference<-if(file.exists(reference_path))read.csv(reference_path,stringsAsFactors=FALSE)else NULL
point_rows<-list();boot_rows<-vector("list",B)

for(m in sizes) {
  ordered<-list(
    Rel=score$miRNA[order(-score$Rel,score$miRNA,method="radix")],
    Pearson=score$miRNA[order(-score$Pearson,score$miRNA,method="radix")],
    RMSE=score$miRNA[order(score$RMSE,score$miRNA,method="radix")])
  sets<-lapply(ordered,function(x)x[seq_len(m)])
  names(sets)<-names(ordered)
  list_path<-file.path(list_dir,paste0("E3_random_lists_m",m,".csv"))
  tab<-read.csv(list_path,stringsAsFactors=FALSE)
  if(!all(c("random_id","miRNA")%in%names(tab))) stop("Malformed random-list file: ",list_path)
  tab<-tab[order(tab$random_id),,drop=FALSE]
  by_id<-split(as.character(tab$miRNA),tab$random_id)
  if(length(by_id)<200L||any(lengths(by_id[seq_len(200L)])!=m)) stop("Need first 200 complete fixed random lists of size ",m)
  random_sets<-lapply(by_id[seq_len(200L)],identity)
  if(any(vapply(c(sets,random_sets),function(x)any(!x%in%colnames(real)),logical(1)))) stop("A selection list contains an unknown feature")
  sets_by_size[[as.character(m)]]<-sets
  random_by_size[[as.character(m)]]<-random_sets
  all_sets<-c(sets,random_sets)
  all_ids<-lapply(all_sets,function(ids)match(ids,colnames(real)))
  rowcors<-lapply(all_ids,function(ix)patient_cor_by_row(real,pred,ix))
  pairs<-network_pairs(all_sets,colnames(real))
  full<-seq_len(nrow(real))
  point_patient<-vapply(rowcors,patient_median,numeric(1),ii=full)
  point_network<-network_values(cor(real,use="pairwise.complete.obs"),cor(pred,use="pairwise.complete.obs"),pairs)
  p_rel<-point_patient[1:3];p_random<-mean(point_patient[4:203],na.rm=TRUE)
  n_rel<-point_network[1:3];n_random<-mean(point_network[4:203],na.rm=TRUE)
  contrasts<-c(Rel_minus_Pearson_patient=unname(p_rel[1]-p_rel[2]),Rel_minus_RMSE_patient=unname(p_rel[1]-p_rel[3]),
    Pearson_minus_Rel_network=unname(n_rel[2]-n_rel[1]),RMSE_minus_Rel_network=unname(n_rel[3]-n_rel[1]),
    Rel_minus_Random_patient=unname(p_rel[1]-p_random),Random_minus_Rel_network=unname(n_random-n_rel[1]))
  point_rows[[as.character(m)]]<-data.frame(size=m,contrast=names(contrasts),estimate=as.numeric(contrasts),stringsAsFactors=FALSE)

  set.seed(20261014L+m)
  for(b in seq_len(B)) {
    ii<-sample.int(nrow(real),nrow(real),replace=TRUE)
    patient_values<-vapply(rowcors,function(rc)patient_median(rc,ii),numeric(1))
    cmr<-suppressWarnings(cor(real[ii,,drop=FALSE],use="pairwise.complete.obs"))
    cmp<-suppressWarnings(cor(pred[ii,,drop=FALSE],use="pairwise.complete.obs"))
    net_values<-network_values(cmr,cmp,pairs)
    pr<-patient_values[1:3];prand<-mean(patient_values[4:203],na.rm=TRUE)
    nr<-net_values[1:3];nrand<-mean(net_values[4:203],na.rm=TRUE)
    vals<-c(Rel_minus_Pearson_patient=unname(pr[1]-pr[2]),Rel_minus_RMSE_patient=unname(pr[1]-pr[3]),
      Pearson_minus_Rel_network=unname(nr[2]-nr[1]),RMSE_minus_Rel_network=unname(nr[3]-nr[1]),
      Rel_minus_Random_patient=unname(pr[1]-prand),Random_minus_Rel_network=unname(nrand-nr[1]))
    if(is.null(boot_rows[[b]])) boot_rows[[b]]<-list()
    boot_rows[[b]][[as.character(m)]]<-vals
  }
  message("Completed paired patient bootstrap: size ",m,"; B=",B,"; seed=",20261014L+m)
}

points<-do.call(rbind,point_rows)
if(!is.null(reference)) {
  ref<-reference[match(paste(points$size,points$contrast),paste(reference$size,reference$contrast)),]
  if(anyNA(ref$estimate)||any(abs(points$estimate-ref$estimate)>1e-10)) stop("Recomputed fixed-set point estimates differ from the existing paired comparison CSV")
}
rows<-list()
for(m in sizes) for(contrast in point_rows[[as.character(m)]]$contrast) {
  vals<-vapply(boot_rows,function(by_size)by_size[[as.character(m)]][[contrast]],numeric(1))
  point<-points$estimate[points$size==m&points$contrast==contrast]
  rows[[length(rows)+1L]]<-data.frame(size=m,contrast=contrast,estimate=point,bootstrap_mean=mean(vals,na.rm=TRUE),
    lower_95=unname(quantile(vals,.025,na.rm=TRUE,type=7)),upper_95=unname(quantile(vals,.975,na.rm=TRUE,type=7)),B=sum(is.finite(vals)))
}
out<-do.call(rbind,rows)
if(B==1000L&&!is.null(reference)) {
  ref<-reference[match(paste(out$size,out$contrast),paste(reference$size,reference$contrast)),]
  if(anyNA(ref$bootstrap_mean)||any(abs(out$bootstrap_mean-ref$bootstrap_mean)>1e-10)||
     any(abs(out$lower_95-ref$lower_95)>1e-10)||any(abs(out$upper_95-ref$upper_95)>1e-10))
    stop("B=1000 bootstrap summary does not reproduce the existing CSV; output was not written")
}
dir.create(dirname(output_path),recursive=TRUE,showWarnings=FALSE)
write.csv(out,output_path,row.names=FALSE)
cat("Wrote ",output_path," (B=",B,")\n",sep="")
