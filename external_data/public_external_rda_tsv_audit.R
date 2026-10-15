load('data/livergabage_cancer_data/chol/data.rda');x<-mrna;y<-mirna
mt<-read.delim('data/livergabage_cancer_data/mRNA.tsv',check.names=FALSE,row.names=1);rt<-read.delim('data/livergabage_cancer_data/mirna.tsv',check.names=FALSE,row.names=1);mp<-read.delim('data/livergabage_cancer_data/anno.probemap',check.names=FALSE)
get_primary_col<-function(pid,cols){z<-cols[substr(cols,1,12)==pid & grepl('-01[A-Z]$',cols)];if(length(z))z[1]else NA_character_}
mi<-intersect(names(y)[-1],rownames(rt)); patient_ids<-x$patient_id
idx<-match(names(x)[-1],mp$gene); ensembl<-mp$id[idx];valid<-!is.na(ensembl)&ensembl%in%rownames(mt)
features<-names(x)[-1][valid];ens<-ensembl[valid]; colmap<-vapply(patient_ids,get_primary_col,character(1),cols=colnames(mt));stopifnot(all(!is.na(colmap)))
a<-as.matrix(x[,features,drop=FALSE]);b<-t(as.matrix(mt[ens,colmap,drop=FALSE]));
cellstats<-function(a,b){a<-as.numeric(a);b<-as.numeric(b);ok<-is.finite(a)&is.finite(b);a<-a[ok];b<-b[ok];tol<-1e-8;data.frame(n=length(a),n_patient=length(patient_ids),n_feature=ncol(as.matrix(x[,features,drop=FALSE])),both_zero=sum(a==0&b==0),both_nonzero=sum(a!=0&b!=0),exact_all=sum(abs(a-b)<tol),exact_both_nonzero=sum(a!=0&b!=0&abs(a-b)<tol),exact_if_either_nonzero=sum((a!=0|b!=0)&abs(a-b)<tol),one_zero=sum((a==0)!=(b==0)),max_abs_difference=max(abs(a-b)),correlation=cor(a,b),rmse=sqrt(mean((a-b)^2)))}
sm<-cellstats(a,b);sm$modality<-'mRNA';sm$transform_match<-'native values against matched Ensembl row from local TSV';
a2<-as.matrix(y[match(patient_ids,y$patient_id),mi,drop=FALSE]);b2<-t(as.matrix(rt[mi,unlist(lapply(patient_ids,get_primary_col,cols=colnames(rt))),drop=FALSE]));sm2<-cellstats(a2,b2);sm2$n_feature<-length(mi);sm2$modality<-'miRNA';sm2$transform_match<-'native values against exact miRNA feature row from local TSV';out<-rbind(sm,sm2);write.csv(out,'new/results/E0/public_external_rda_vs_local_tsv.csv',row.names=FALSE);print(out)
