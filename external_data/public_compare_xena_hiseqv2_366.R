load('data/data.rda')
stopifnot(exists('mrna'), nrow(mrna)==366L)
xena_path <- 'new/results/E0/public_TCGA_LIHC_HiSeqV2.gz'
xena <- read.delim(gzfile(xena_path), header=TRUE, row.names=1, check.names=FALSE,
                    quote='', comment.char='', stringsAsFactors=FALSE)
patients <- unique(as.character(mrna$patient_id))
features <- setdiff(names(mrna), 'patient_id')
shared <- intersect(features, rownames(xena))
if(length(patients)!=366L) stop('Expected 366 unique patient IDs')
if(length(shared)!=length(features)) warning('Not all input genes are in Xena')
map <- data.frame(patient_id=patients, expected_sample=paste0(patients,'-01'),
                  candidate_count=0L, selected_sample=NA_character_, status='missing',
                  stringsAsFactors=FALSE)
summary <- vector('list', length(patients)); diffs <- list()
for(i in seq_along(patients)) {
  pid <- patients[i]; target <- paste0(pid,'-01')
  exact_cols <- which(colnames(xena)==target)
  map$candidate_count[i] <- length(exact_cols)
  if(length(exact_cols)!=1L) {
    map$status[i] <- if(length(exact_cols)==0L) 'missing' else 'duplicate_exact_sample'
    summary[[i]] <- data.frame(patient_id=pid, selected_sample=NA_character_, n_shared_genes=length(shared),
      n_compared=0L,n_exact=NA_integer_,n_within_1e_6=NA_integer_,n_nonzero_both=NA_integer_,
      n_nonzero_both_within_1e_6=NA_integer_,max_abs_diff=NA_real_,rmse=NA_real_,pearson=NA_real_)
    next
  }
  j <- exact_cols[1]; map$selected_sample[i] <- colnames(xena)[j]; map$status[i] <- 'matched'
  row <- mrna[mrna$patient_id==pid,,drop=FALSE]
  if(nrow(row)!=1L) stop('Non-unique RDA case: ',pid)
  a <- as.numeric(row[1,shared]); b <- as.numeric(xena[shared,j]); ok <- is.finite(a)&is.finite(b)
  d <- b-a; nz <- ok&a!=0&b!=0
  summary[[i]] <- data.frame(patient_id=pid,selected_sample=colnames(xena)[j],n_shared_genes=length(shared),
     n_compared=sum(ok),n_exact=sum(ok&d==0),n_within_1e_6=sum(ok&abs(d)<1e-6),
     n_nonzero_both=sum(nz),n_nonzero_both_within_1e_6=sum(nz&abs(d)<1e-6),
     max_abs_diff=max(abs(d[ok])),rmse=sqrt(mean(d[ok]^2)),pearson=cor(a[ok],b[ok]))
  bad <- which(ok & abs(d)>=1e-6)
  if(length(bad)) diffs[[length(diffs)+1L]] <- data.frame(patient_id=pid,sample_id=colnames(xena)[j],gene=shared[bad],
      input_value=a[bad],xena_value=b[bad],difference_xena_minus_input=d[bad],stringsAsFactors=FALSE)
}
summary <- do.call(rbind,summary)
write.csv(map,'new/results/E0/public_LIHC_Xena_HiSeqV2_case_map_366.csv',row.names=FALSE,na='')
write.csv(summary,'new/results/E0/public_LIHC_Xena_HiSeqV2_match_summary_366.csv',row.names=FALSE,na='')
if(length(diffs)) write.csv(do.call(rbind,diffs),'new/results/E0/public_LIHC_Xena_HiSeqV2_differences_366.csv',row.names=FALSE) else write.csv(data.frame(patient_id=character(),sample_id=character(),gene=character(),input_value=double(),xena_value=double(),difference_xena_minus_input=double()),'new/results/E0/public_LIHC_Xena_HiSeqV2_differences_366.csv',row.names=FALSE)
cat('patients=',length(patients),'matched=',sum(map$status=='matched'),'missing=',sum(map$status=='missing'),'duplicate_exact=',sum(map$status=='duplicate_exact_sample'),'shared_genes=',length(shared),'rows_diff=',sum(summary$n_within_1e_6<summary$n_compared,na.rm=TRUE),'\n')
print(summary(summary$n_within_1e_6/summary$n_compared))
