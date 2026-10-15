load('data/data.rda'); d<-mrna; d<-d[d$patient_id=='TCGA-XR-A8TC',,drop=FALSE]
g<-read.delim('new/results/E0/public_gdc_sample_data/lihc_XR-A8TC_star.tsv',check.names=FALSE,comment.char='#');g<-g[!is.na(g$gene_name)&nzchar(g$gene_name),];g<-g[!duplicated(g$gene_name),]
features<-names(d)[-1];j<-match(features,g$gene_name);ok<-!is.na(j); y<-as.numeric(d[1,features[ok]])
transforms<-list(unstranded_counts=as.numeric(g$unstranded[j[ok]]),TPM=as.numeric(g$tpm_unstranded[j[ok]]),FPKM=as.numeric(g$fpkm_unstranded[j[ok]]),FPKM_UQ=as.numeric(g$fpkm_uq_unstranded[j[ok]]))
out<-list();for(nm in names(transforms)){x<-transforms[[nm]];for(tr in c('identity','log2_1p','ln_1p')){z<-if(tr=='identity')x else if(tr=='log2_1p')log2(x+1)else log(x+1); good<-is.finite(z)&is.finite(y);out[[length(out)+1]]<-data.frame(source=nm,transformation=tr,n=sum(good),cor=cor(z[good],y[good]),rmse=sqrt(mean((z[good]-y[good])^2)),exact_fraction=mean(abs(z[good]-y[good])<1e-6))}}
write.csv(do.call(rbind,out),'new/results/E0/public_LIHC_mRNA_transform_checks.csv',row.names=FALSE);print(do.call(rbind,out))
