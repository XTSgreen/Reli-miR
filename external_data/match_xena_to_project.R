# Match each Xena matrix's patient set to a GDC project.
basics <- read.csv("new/data/xena_file_basics.csv", stringsAsFactors = FALSE)
projects <- c("LIHC", "KIRC", "LUAD", "COAD", "STAD", "PRAD", "THCA", "HNSC")
gdc <- setNames(lapply(projects, function(p)
  readLines(sprintf("new/data/gdc_cases_TCGA_%s.txt", p))), projects)
rows <- lapply(seq_len(nrow(basics)), function(i) {
  pf <- sprintf("new/data/patients_%s.txt",
    gsub("[^A-Za-z0-9]", "_", basics$file[i]))
  pats <- readLines(pf)
  hits <- sapply(gdc, function(g) length(intersect(pats, g)))
  out <- data.frame(file = basics$file[i], n_patients = length(pats),
    stringsAsFactors = FALSE)
  for (p in projects) out[[p]] <- as.integer(hits[[p]])
  out$assigned_project <- projects[which.max(hits)]
  out$exclusive <- sum(hits > 0) == 1L
  out
})
out <- do.call(rbind, rows)
write.csv(out, "new/data/xena_file_cohort_map.csv", row.names = FALSE)
print(out)
