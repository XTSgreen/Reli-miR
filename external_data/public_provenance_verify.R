# Reproducible, local-only value checks for the supplied RDA matrices and candidate TCGA-CHOL public files.
load('data/data.rda')
primary <- clinical
write.csv(primary, 'new/results/E0/public_primary_clinical.csv', row.names=FALSE, na='')
load('data/livergabage_cancer_data/chol/data.rda')
external_mrna <- mrna
external_mirna <- mirna
# Read local Xena/GDC-style candidate files without uploading any data.
mrna_tsv <- read.delim('data/livergabage_cancer_data/mRNA.tsv', check.names=FALSE, row.names=1)
mirna_tsv <- read.delim('data/livergabage_cancer_data/mirna.tsv', check.names=FALSE, row.names=1)
map <- read.delim('data/livergabage_cancer_data/anno.probemap', check.names=FALSE)
# Match samples and features for the gene matrix via the local probemap.
primary_for_patient <- function(ids, sample_cols) {
  out <- vapply(ids, function(pid) {
    hits <- sample_cols[substr(sample_cols, 1, 12) == pid & grepl('-01[A-Z]$', sample_cols)]
    if (length(hits)) hits[[1]] else NA_character_
  }, character(1))
  names(out) <- ids
  out
}
sample_for_m <- primary_for_patient(external_mrna$patient_id, colnames(mrna_tsv))
sample_for_mi <- primary_for_patient(external_mirna$patient_id, colnames(mirna_tsv))
common_samples_m <- names(sample_for_m)[!is.na(sample_for_m)]
external_gene <- names(external_mrna)[-1]
map_gene <- map[match(external_gene, map$gene), , drop=FALSE]
common_gene_rows <- which(!is.na(map_gene$id) & map_gene$id %in% rownames(mrna_tsv))
# Match miRNA rows by exact precursor-style identifier.
common_samples_mi <- names(sample_for_mi)[!is.na(sample_for_mi)]
common_mi <- intersect(names(external_mirna)[-1], rownames(mirna_tsv))
cor_compare <- function(a,b,transform, label) {
  if (transform=='identity') z <- a else if(transform=='log2_1p') z <- log2(pmax(a,0)+1) else if(transform=='ln_1p') z <- log(pmax(a,0)+1)
  data.frame(transform=transform,n=length(a),exact_fraction=mean(abs(z-b)<1e-8,na.rm=TRUE),max_abs_diff=max(abs(z-b),na.rm=TRUE),correlation=cor(z,b,use='complete.obs'))
}
# Use every exact shared sample-feature cell for miRNA; gene rows map Ensembl IDs to symbols.
rows <- list()
if(length(common_samples_mi)>0 && length(common_mi)>0) {
  a <- as.matrix(external_mirna[match(common_samples_mi,external_mirna$patient_id),common_mi,drop=FALSE])
  b <- t(as.matrix(mirna_tsv[common_mi,unname(sample_for_mi[common_samples_mi]),drop=FALSE]))
  rows[[length(rows)+1]] <- transform(cor_compare(as.numeric(a),as.numeric(b),'identity','miRNA'))
}
if(length(common_samples_m)>0 && length(common_gene_rows)>0) {
  genes <- external_gene[common_gene_rows]
  ensembl <- map_gene$id[common_gene_rows]
  a <- as.matrix(external_mrna[match(common_samples_m,external_mrna$patient_id),genes,drop=FALSE])
  b <- t(as.matrix(mrna_tsv[ensembl,unname(sample_for_m[common_samples_m]),drop=FALSE]))
  for (tr in c('identity','log2_1p','ln_1p')) rows[[length(rows)+1]] <- cor_compare(as.numeric(a),as.numeric(b),tr,'mRNA')
}
if(length(rows)) write.csv(do.call(rbind,rows),'new/results/E0/public_external_expression_match.csv',row.names=FALSE)
# Per-sample survival comparison with the adjacent clinical and survival tables.
clin_tsv <- read.delim('data/livergabage_cancer_data/clinical.tsv',check.names=FALSE)
surv_tsv <- read.delim('data/livergabage_cancer_data/survival.tsv',check.names=FALSE)
get_patient <- function(x) substr(x,1,12)
# Convert patient-level metadata from GDC/Xena clinical file, preserving unavailable follow-up values.
clin_tsv$patient <- get_patient(clin_tsv$sample)
surv_tsv$patient <- get_patient(surv_tsv$sample)
# Save candidate source field values; exact primary OS comparison uses GDC API output in the companion verifier.
write.csv(unique(clin_tsv[clin_tsv$project_id.project=='TCGA-CHOL', c('patient','project_id.project','vital_status.demographic','days_to_death.demographic','days_to_last_follow_up.diagnoses')]), 'new/results/E0/public_external_clinical_fields.csv', row.names=FALSE, na='')
write.csv(surv_tsv, 'new/results/E0/public_external_survival_rows.csv', row.names=FALSE, na='')
cat('mRNA matched samples',length(common_samples_m),'mapped genes',length(common_gene_rows),'\n')
cat('miRNA matched samples',length(common_samples_mi),'features',length(common_mi),'\n')
cat('Wrote public-source candidate comparisons.\n')
