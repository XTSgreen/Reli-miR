source_dir <- 'new/results/E4_harmonized_inputs/source_cache'
out_dir <- 'new/results/E4_harmonized_inputs'
input_rda <- 'data/livergabage_cancer_data/chol/data.rda'
source_xena <- file.path(source_dir,'public_TCGA_CHOL_HiSeqV2.gz')
source_meta <- file.path(source_dir,'public_TCGA_CHOL_HiSeqV2.json')
stopifnot(file.exists(source_xena),file.exists(source_meta))
orig <- new.env(parent=emptyenv()); load(input_rda,envir=orig)
stopifnot(all(c('mrna','mirna','clinical') %in% ls(orig)))
patients <- as.character(orig$mrna$patient_id)
stopifnot(length(patients)==35L, !anyDuplicated(patients))
xena <- read.delim(gzfile(source_xena),header=TRUE,row.names=1,check.names=FALSE,
                   quote='',comment.char='',stringsAsFactors=FALSE)
if(anyDuplicated(rownames(xena))) stop('Duplicate Xena feature IDs; refusing ambiguous mapping')
targets <- paste0(patients,'-01')
candidate_count <- vapply(targets,function(s)sum(colnames(xena)==s),integer(1))
map <- data.frame(patient_id=patients,expected_primary_sample=targets,
                  exact_candidate_count=candidate_count,
                  selected_sample=ifelse(candidate_count==1L,targets,NA_character_),
                  status=ifelse(candidate_count==1L,'matched',ifelse(candidate_count==0L,'missing','ambiguous_duplicate')),
                  stringsAsFactors=FALSE)
write.csv(map,file.path(out_dir,'CHOL_35_case_sample_mapping.csv'),row.names=FALSE,na='')
if(any(candidate_count!=1L)) stop('Every case must map to exactly one -01 Xena sample; see mapping CSV')
model_genes <- read.csv('new/results/E0/E0_model_input_gene_list.csv',stringsAsFactors=FALSE,check.names=FALSE)$gene
if(length(model_genes)!=3000L || anyDuplicated(model_genes)) stop('Unexpected model gene list')
xena_genes <- rownames(xena)
coverage <- data.frame(gene=model_genes,in_xena=model_genes %in% xena_genes,
                       stringsAsFactors=FALSE)
write.csv(coverage,file.path(out_dir,'CHOL_E0_3000_gene_coverage.csv'),row.names=FALSE)
selected <- xena[,targets,drop=FALSE]
mrna <- as.data.frame(t(as.matrix(selected)),check.names=FALSE,stringsAsFactors=FALSE)
mrna <- data.frame(patient_id=patients,mrna,check.names=FALSE,stringsAsFactors=FALSE)
mirna <- orig$mirna
clinical <- orig$clinical
stopifnot(nrow(mrna)==35L,ncol(mrna)==nrow(xena)+1L,
          setequal(as.character(mirna$patient_id),patients),
          nrow(clinical)==0L,ncol(clinical)==85L)
out_rda <- file.path(out_dir,'external_CHOL_matched_RSEM.rda')
save(mrna,mirna,clinical,file=out_rda,version=2)
cat('cases=',nrow(mrna),'features=',ncol(mrna)-1L,'model_genes=',length(model_genes),
    'covered=',sum(coverage$in_xena),'missing=',sum(!coverage$in_xena),
    'mirna=',nrow(mirna),'x',ncol(mirna)-1L,'clinical=',nrow(clinical),'x',ncol(clinical),'\n')

