# Identify which TCGA cohort each downloaded Xena HiSeqV2 matrix belongs to,
# by matching patient barcodes against GDC case lists per project.
# Writes: new/data/xena_file_basics.csv and new/data/patients_<file>.txt
files <- c(
  "new/data/HiSeqV2.gz",
  "new/data/HiSeqV2 (1).gz",
  "new/data/HiSeqV2 (2).gz",
  "new/data/HiSeqV2 (3).gz",
  "new/data/HiSeqV2 (4).gz",
  "new/data/HiSeqV2 (5).gz",
  "new/data/HiSeqV2 (6).gz")

info <- do.call(rbind, lapply(files, function(f) {
  con <- gzfile(f, "rt"); on.exit(close(con))
  hdr <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
  rest <- readLines(con)
  patients <- unique(substr(tail(hdr, -1), 1, 12))
  writeLines(patients, sprintf("new/data/patients_%s.txt",
    gsub("[^A-Za-z0-9]", "_", basename(f))))
  data.frame(file = basename(f), n_cols = length(hdr), n_genes = length(rest),
             n_patients = length(patients), stringsAsFactors = FALSE)
}))
write.csv(info, "new/data/xena_file_basics.csv", row.names = FALSE)
print(info)
