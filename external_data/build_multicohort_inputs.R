# Build same-pipeline external inputs for 7 TCGA cohorts.
# mRNA: Xena HiSeqV2 log2(RSEM norm_count+1) (primary-tumor -01 samples).
# miRNA: GDC mirbase21 quantification, RPM -> log2(RPM+1), hairpin ID space.
# Output: new/results/E4_multicohort/inputs/external_<C>_matched_RSEM.rda
# plus validation CSVs (sample mapping, gene coverage, miRNA overlap).
suppressWarnings(suppressMessages(library(readr)))
cohorts <- c("KIRC", "LUAD", "COAD", "STAD", "PRAD", "THCA", "HNSC")
cohort_map <- read.csv("new/data/xena_file_cohort_map.csv", stringsAsFactors = FALSE)
canon <- readRDS("new/results/E0/canonical_input.rds")
internal_mirnas <- canon$miRNA_names
stopifnot(!is.null(internal_mirnas), length(internal_mirnas) > 0)
model_genes <- read.csv("new/results/E0/E0_model_input_gene_list.csv",
                        stringsAsFactors = FALSE, check.names = FALSE)$gene
stopifnot(length(model_genes) == 3000L)
chol_env <- new.env(parent = emptyenv())
load("data/livergabage_cancer_data/chol/data.rda", envir = chol_env)
clin_names <- names(chol_env$clinical)

read_quant <- function(p) {
  magic <- readBin(p, "raw", n = 2L)
  con <- if (length(magic) == 2L && magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b))
    gzfile(p, "rt") else file(p, "rt")
  on.exit(close(con))
  tab <- read.delim(con, header = TRUE, check.names = FALSE,
                    quote = "", comment.char = "", stringsAsFactors = FALSE)
  tab <- tab[, c("miRNA_ID", "reads_per_million_miRNA_mapped"), drop = FALSE]
  names(tab) <- c("miRNA", "RPM")
  tab
}

read_xena <- function(p) {
  x <- read.delim(gzfile(p), header = TRUE, row.names = 1, check.names = FALSE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
  x
}

summary_rows <- list()
for (co in cohorts) {
  out_dir <- file.path("new/results/E4_multicohort/inputs", co)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(out_dir, sprintf("external_%s_matched_RSEM.rda", co)))) {
    cat(co, "input already built - skipping\n"); next
  }
  man <- read.csv(sprintf("new/data/gdc_mirna/manifest_%s.csv", co),
                  stringsAsFactors = FALSE, colClasses = "character")
  man$path <- file.path("new/data/gdc_mirna", paste0("raw_", co), man$file_id, man$file_name)
  man$exists <- file.exists(man$path)
  if (sum(man$exists) < 0.99 * nrow(man)) {
    cat(co, "raw files incomplete (", sum(man$exists), "/", nrow(man),
        "); skipping until download finishes\n")
    next
  }
  man <- man[order(man$case_id, man$sample_id, man$file_id), ]
  keep <- !duplicated(man$case_id)
  map <- transform(man, status = ifelse(!exists, "file_missing",
    ifelse(keep, "matched", "dropped_duplicate_case_file")))
  sel <- man[keep & man$exists, , drop = FALSE]
  write.csv(map[, c("case_id", "sample_id", "file_id", "status")],
            file.path(out_dir, "sample_mapping.csv"), row.names = FALSE, na = "")
  cat(co, "manifest", nrow(man), "existing", sum(man$exists),
      "selected_cases", nrow(sel), "\n")

  mirna_list <- lapply(seq_len(nrow(sel)), function(i) {
    q <- read_quant(sel$path[i])
    setNames(log2(pmax(q$RPM, 0) + 1), q$miRNA)
  })
  names(mirna_list) <- sel$case_id
  mirna_mat <- do.call(rbind, lapply(mirna_list, function(v) as.data.frame(t(v))))
  mirna_mat[is.na(mirna_mat)] <- NA_real_

  mrna_file <- if (co == "COAD") "new/data/xena_mrna_TCGA_COAD_HiSeqV2.gz" else
    file.path("new/data", cohort_map$file[match(co, cohort_map$assigned_project)])
  xz <- read_xena(mrna_file)
  cols01 <- grep("-01$", colnames(xz), value = TRUE)
  stopifnot(length(cols01) > 0)
  pats_xena <- substr(cols01, 1, 12)
  if (anyDuplicated(pats_xena)) {
    dup <- unique(pats_xena[duplicated(pats_xena)])
    keep_col <- !(pats_xena %in% dup)
    pats_xena <- pats_xena[keep_col]; cols01 <- cols01[keep_col]
  } else keep_col <- rep(TRUE, length(cols01))
  common_pat <- intersect(sel$case_id, pats_xena)
  stopifnot(length(common_pat) > 20)
  mi_col <- match(common_pat, sel$case_id)
  mx_col <- match(common_pat, pats_xena)
  mrna_sub <- t(as.matrix(xz[, cols01[mx_col], drop = FALSE]))
  rownames(mrna_sub) <- common_pat
  colnames(mrna_sub) <- rownames(xz)
  mirna_sub <- mirna_mat[common_pat, , drop = FALSE]
  cov <- model_genes %in% colnames(mrna_sub)
  write.csv(data.frame(gene = model_genes, in_external = cov),
            file.path(out_dir, "gene_coverage_3000.csv"), row.names = FALSE)
  mirna_overlap <- intersect(internal_mirnas, colnames(mirna_sub))
  write.csv(data.frame(internal_miRNA = internal_mirnas,
                       in_external = internal_mirnas %in% colnames(mirna_sub)),
            file.path(out_dir, "miRNA_overlap_846.csv"), row.names = FALSE)
  mrna_df <- data.frame(patient_id = rownames(mrna_sub), mrna_sub,
                        check.names = FALSE, stringsAsFactors = FALSE)
  mirna_df <- data.frame(patient_id = rownames(mirna_sub), mirna_sub,
                         check.names = FALSE, stringsAsFactors = FALSE)
  clinical <- as.data.frame(matrix(logical(0), nrow = 0, ncol = length(clin_names),
                                   dimnames = list(NULL, clin_names)))
  stopifnot(nrow(mrna_df) == length(common_pat), nrow(mirna_df) == length(common_pat))
  mrna <- mrna_df; mirna <- mirna_df
  save(mrna, mirna, clinical,
       file = file.path(out_dir, sprintf("external_%s_matched_RSEM.rda", co)), version = 2)
  summary_rows[[co]] <- data.frame(cohort = co, n_patients = length(common_pat),
    n_genes_xena = ncol(mrna_sub), gene_coverage_3000 = mean(cov),
    n_mirna_external = ncol(mirna_sub), n_overlap_internal_846 = length(mirna_overlap),
    mirna_value_min = min(mirna_sub, na.rm = TRUE), mirna_value_max = max(mirna_sub, na.rm = TRUE),
    mrna_value_min = min(mrna_sub, na.rm = TRUE), mrna_value_max = max(mrna_sub, na.rm = TRUE))
  cat(co, "DONE patients", length(common_pat), "genes covered", sum(cov), "/3000; mirna cols",
      ncol(mirna_sub), "overlap", length(mirna_overlap), "\n")
}
sum_df <- do.call(rbind, summary_rows)
write.csv(sum_df, "new/results/E4_multicohort/inputs/input_build_summary.csv", row.names = FALSE)
print(sum_df)
