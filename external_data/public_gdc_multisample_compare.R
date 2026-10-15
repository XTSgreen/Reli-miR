# Reproducible, local-only validation of supplied matrices against per-case files from the official GDC API.
load('data/data.rda'); primary_mrna<-mrna; primary_miRNA<-mirna
load('data/livergabage_cancer_data/chol/data.rda'); chol_mrna<-mrna; chol_miRNA<-mirna
samples<-data.frame(cohort=c('CHOL','CHOL','CHOL','LIHC','LIHC'),patient=c('TCGA-W5-AA2X','TCGA-W5-AA30','TCGA-ZH-A8Y2','TCGA-XR-A8TC','TCGA-BC-A112'),star=c('chol_01A_star_counts.tsv','chol_W5-AA30_star.tsv','chol_ZH-A8Y2_star.tsv','lihc_XR-A8TC_star.tsv','lihc_BC-A112_star.tsv'),mirna=c('chol_01A_mirna_quantification.txt','chol_W5-AA30_miRNA.txt','chol_ZH-A8Y2_miRNA.txt','lihc_XR-A8TC_miRNA.txt','lihc_BC-A112_miRNA.txt'),stringsAsFactors=FALSE)
compare<-function(x,y){ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok];data.frame(n=length(x),correlation=cor(x,y),rmse=sqrt(mean((x-y)^2)),fraction_within_1e_6=mean(abs(x-y)<1e-6),max_abs_difference=max(abs(x-y)))}
out<-list()
for(i in seq_len(nrow(samples))){
 s<-samples[i,]; is_chol<-s$cohort=='CHOL'; xm<-if(is_chol) chol_mrna else primary_mrna; xr<-if(is_chol) chol_miRNA else primary_miRNA
 star<-read.delim(file.path('new/results/E0/public_gdc_sample_data',s$star),check.names=FALSE,comment.char='#')
 mir<-read.delim(file.path('new/results/E0/public_gdc_sample_data',s$mirna),check.names=FALSE)
 # GDC per-case miRNA file reports raw read counts and RPM; the candidate conversion tested is log2(RPM+1).
 ridx<-match(s$patient,xr$patient_id); mir_features<-names(xr)[-1]; mj<-match(mir_features,mir$miRNA_ID); ok<-!is.na(mj)
 observed<-as.numeric(xr[ridx,mir_features[ok]]); expected<-log2(as.numeric(mir$reads_per_million_miRNA_mapped[mj[ok]])+1)
 q<-compare(expected,observed); out[[length(out)+1]]<-data.frame(cohort=s$cohort,patient=s$patient,modality='miRNA',tested_transform='log2(RPM+1)',q)
 # GDC v36 symbol mapping is not one-to-one; use the first row per gene_name without choosing by value.
 star<-star[!is.na(star$gene_name)&nzchar(star$gene_name),]; star<-star[!duplicated(star$gene_name),]
 features<-names(xm)[-1]; gi<-match(features,star$gene_name); ok<-!is.na(gi)
 observed<-as.numeric(xm[match(s$patient,xm$patient_id),features[ok]])
 expected<-log2(as.numeric(star$unstranded[gi[ok]])+1)
 q<-compare(expected,observed); out[[length(out)+1]]<-data.frame(cohort=s$cohort,patient=s$patient,modality='mRNA',tested_transform='log2(unstranded STAR counts+1)',q)
}
res<-do.call(rbind,out);write.csv(res,'new/results/E0/public_gdc_multisample_expression_match.csv',row.names=FALSE)
print(res)
