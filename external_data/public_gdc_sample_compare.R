load('data/livergabage_cancer_data/chol/data.rda'); sid <- 'TCGA-W5-AA2X'; mi <- mrna[mrna$patient_id==sid,,drop=FALSE]; mirna_one <- mirna[mirna$patient_id==sid,,drop=FALSE]
star <- read.delim('new/results/E0/public_gdc_sample_data/chol_01A_star_counts.tsv',check.names=FALSE,comment.char='#')
mir <- read.delim('new/results/E0/public_gdc_sample_data/chol_01A_mirna_quantification.txt',check.names=FALSE)
compare <- function(x,y) c(correlation=cor(x,y,use='pairwise.complete.obs'),RMSE=sqrt(mean((x-y)^2,na.rm=TRUE)),fraction_within_1e_6=mean(abs(x-y)<1e-6,na.rm=TRUE))
out <- list()
# Match by exact legacy GDC miRNA precursor ID; compare current GDC read counts and RPM.
mi_ids <- names(mirna_one)[-1]
j <- match(mi_ids,mir$miRNA_ID); ok <- !is.na(j)
obs <- as.numeric(mirna_one[1,mi_ids[ok]])
for(key in c('read_count','reads_per_million_miRNA_mapped')) {
 x<-as.numeric(mir[[key]][j[ok]])
 for(tr in c('identity','log2_1p','ln_1p')) {
  z<-if(tr=='identity') x else if(tr=='log2_1p') log2(x+1) else log(x+1)
  q<-compare(z,obs); out[[length(out)+1]]<-data.frame(modality='miRNA',source=key,transformation=tr,n=sum(is.finite(z)&is.finite(obs)),correlation=q[1],rmse=q[2],fraction_within_1e_6=q[3])
 }
}
# Match symbols to the first GDC gene row per gene name; this is a diagnostic only because gene naming versions differ.
gene_tab <- star[!is.na(star$gene_name)&nzchar(star$gene_name),]; gene_tab<-gene_tab[!duplicated(gene_tab$gene_name),]
gene_names <- names(mi)[-1]; j<-match(gene_names,gene_tab$gene_name); ok<-!is.na(j)
obs<-as.numeric(mi[1,gene_names[ok]]); gdc<-gene_tab[j[ok],]
for(key in c('unstranded','tpm_unstranded','fpkm_unstranded','fpkm_uq_unstranded')) {
 x<-as.numeric(gdc[[key]])
 for(tr in c('identity','log2_1p','ln_1p')) {
  z<-if(tr=='identity') x else if(tr=='log2_1p') log2(x+1) else log(x+1)
  q<-compare(z,obs); out[[length(out)+1]]<-data.frame(modality='mRNA',source=key,transformation=tr,n=sum(is.finite(z)&is.finite(obs)),correlation=q[1],rmse=q[2],fraction_within_1e_6=q[3])
 }
}
write.csv(do.call(rbind,out),'new/results/E0/public_gdc_sample_value_match.csv',row.names=FALSE)
cat('miRNA matched',sum(ok),'; full comparison table written.\n')

