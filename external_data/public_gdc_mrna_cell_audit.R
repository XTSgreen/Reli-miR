suppressPackageStartupMessages(library(jsonlite))
read_j<-function(path){s<-rawToChar(readBin(path,'raw',n=file.info(path)$size));fromJSON(sub('^\\ufeff','',s),simplifyVector=FALSE)}
getv<-function(x)if(is.null(x)||!length(x))numeric()else suppressWarnings(as.numeric(unlist(x,use.names=FALSE)))
collect<-function(path)read_j(path)$data$hits
compare_cells<-function(x,y){ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok];tol<-1e-6;data.frame(n=length(x),both_zero=sum(x==0&y==0),both_nonzero=sum(x!=0&y!=0),exact_all=sum(abs(x-y)<tol),exact_both_nonzero=sum(x!=0&y!=0&abs(x-y)<tol),exact_if_either_nonzero=sum((x!=0|y!=0)&abs(x-y)<tol),nonzero_disagree=sum((x==0)!=(y==0)),rmse=sqrt(mean((x-y)^2)),cor=cor(x,y))}
load('data/data.rda'); pm<-mrna;pmi<-mirna;primary<-clinical
load('data/livergabage_cancer_data/chol/data.rda');cm<-mrna;cmi<-mirna
patients<-list(CHOL=cm$patient_id,LIHC=c('TCGA-XR-A8TC','TCGA-BC-A112'))
filemap<-data.frame(patient=c('TCGA-W5-AA2X','TCGA-W5-AA30','TCGA-ZH-A8Y2','TCGA-XR-A8TC','TCGA-BC-A112'),star=c('chol_01A_star_counts.tsv','chol_W5-AA30_star.tsv','chol_ZH-A8Y2_star.tsv','lihc_XR-A8TC_star.tsv','lihc_BC-A112_star.tsv'),mir=c('chol_01A_mirna_quantification.txt','chol_W5-AA30_miRNA.txt','chol_ZH-A8Y2_miRNA.txt','lihc_XR-A8TC_miRNA.txt','lihc_BC-A112_miRNA.txt'))
# Detail only mRNA; count transformed public official current STAR field compared with source native matrix.
out<-list()
for(i in 1:nrow(filemap)){
 pid<-filemap$patient[i]; chol<-i<=3;dat<-if(chol)cm else pm
 g<-read.delim(file.path('new/results/E0/public_gdc_sample_data',filemap$star[i]),check.names=FALSE,comment.char='#');g<-g[!is.na(g$gene_name)&nzchar(g$gene_name),];g<-g[!duplicated(g$gene_name),]
 r<-match(pid,dat$patient_id); names0<-names(dat)[-1];j<-match(names0,g$gene_name);ok<-!is.na(j); obs<-as.numeric(dat[r,names0[ok]]); raw<-as.numeric(g$unstranded[j[ok]]); exp<-log2(raw+1); z<-compare_cells(obs,exp); z$cohort<-if(chol)'CHOL'else'LIHC';z$patient<-pid; z$tested_gene_symbols<-sum(ok);out[[i]]<-z
}
res<-do.call(rbind,out);write.csv(res,'new/results/E0/public_gdc_mrna_match_cell_audit.csv',row.names=FALSE);print(res)
# Ensembl gene name duplication counts among the gene symbols used for the RDA.
