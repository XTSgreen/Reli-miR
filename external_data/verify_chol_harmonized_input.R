orig <- new.env(parent=emptyenv()); load('data/livergabage_cancer_data/chol/data.rda',envir=orig)
harm <- new.env(parent=emptyenv()); load('new/results/E4_harmonized_inputs/external_CHOL_matched_RSEM.rda',envir=harm)
xena <- read.delim(gzfile('new/results/E4_harmonized_inputs/source_cache/public_TCGA_CHOL_HiSeqV2.gz'),header=TRUE,row.names=1,check.names=FALSE,quote='',comment.char='',stringsAsFactors=FALSE)
patients <- as.character(orig$mrna$patient_id); features <- rownames(xena); targets <- paste0(patients,'-01')
stopifnot(nrow(orig$mrna)==35L, all(vapply(targets,function(x)sum(colnames(xena)==x),integer(1))==1L))
stopifnot(identical(harm$mirna,orig$mirna),identical(harm$clinical,orig$clinical))
stopifnot(identical(as.character(harm$mrna$patient_id),patients),identical(names(harm$mrna)[-1],features))
maxdiff <- 0; ncomp <- 0L
for(i in seq_along(patients)) {
  a <- as.numeric(harm$mrna[i,-1]); b <- as.numeric(xena[,targets[i]])
  stopifnot(length(a)==length(b),all(is.finite(a)),all(is.finite(b)))
  maxdiff <- max(maxdiff,abs(a-b)); ncomp <- ncomp+length(a)
}
stopifnot(maxdiff==0,ncomp==35L*length(features))
model_genes <- read.csv('new/results/E0/E0_model_input_gene_list.csv',stringsAsFactors=FALSE)$gene
coverage <- data.frame(model_gene_count=length(model_genes),covered=sum(model_genes%in%features),missing=sum(!model_genes%in%features))
write.csv(data.frame(input_matrix_values_compared=ncomp,max_abs_difference=maxdiff,
                     cases=nrow(harm$mrna),xena_features=length(features),mirna_unchanged=TRUE,
                     clinical_unchanged=TRUE,model_gene_count=length(model_genes),
                     model_genes_covered=sum(model_genes%in%features),model_genes_missing=sum(!model_genes%in%features)),
          'new/results/E4_harmonized_inputs/CHOL_harmonized_input_validation.csv',row.names=FALSE)
cat('validated cases=',nrow(harm$mrna),'values=',ncomp,'max_abs_diff=',maxdiff,
    'mirna_identical=TRUE clinical_identical=TRUE model_coverage=',sum(model_genes%in%features),'/',length(model_genes),'\n')
