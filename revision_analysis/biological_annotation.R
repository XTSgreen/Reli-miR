#!/usr/bin/env Rscript
# E5 descriptive annotation and expression analysis for canonical precursor miRNA loci.
# Run from repository root: Rscript new/code/biological_annotation.R
options(stringsAsFactors = FALSE)
out <- "new/results/E5"
ann <- file.path(out, "annotation")
dir.create(ann, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
can <- readRDS("new/results/E0/canonical_input.rds")
scores <- read.csv("new/results/E0/PCA_ElasticNet_reliability_from_cal.csv", check.names = FALSE)
stopifnot(nrow(scores) == length(can$miRNA_names), !anyDuplicated(scores$miRNA),
          setequal(scores$miRNA, can$miRNA_names))
raw_env <- new.env(parent = emptyenv())
load("data/data.rda", envir = raw_env)
raw_mrna <- raw_env$mrna
all_y_native <- rbind(can$Y_train_native, can$Y_cal_native, can$Y_test_native)
all_patient_ids <- rownames(all_y_native)
raw_gene_names <- setdiff(names(raw_mrna), "patient_id")
stopifnot(nrow(raw_mrna) == 366L, length(raw_gene_names) == 18204L,
          !anyDuplicated(raw_mrna$patient_id), !anyDuplicated(raw_gene_names),
          setequal(as.character(raw_mrna$patient_id), all_patient_ids),
          all(can$gene_names %in% raw_gene_names))
all_mrna <- as.matrix(raw_mrna[match(all_patient_ids, raw_mrna$patient_id), raw_gene_names, drop = FALSE])
storage.mode(all_mrna) <- "double"
rownames(all_mrna) <- all_patient_ids
model_input_values <- all_mrna[, can$gene_names, drop = FALSE]
canonical_model_values <- rbind(can$X_train_native, can$X_cal_native, can$X_test_native)
model_value_delta <- max(abs(model_input_values - canonical_model_values), na.rm = TRUE)
stopifnot(identical(rownames(model_input_values), rownames(canonical_model_values)),
          model_value_delta < 1e-10)
raw_mrna_alignment <- list(n_patients = nrow(all_mrna), n_mrna_features = ncol(all_mrna),
  n_model_input_features = ncol(model_input_values), model_feature_max_abs_difference_vs_canonical = model_value_delta,
  alignment = "raw mRNA ordered by exact patient_id match to concatenated canonical train/cal/test Y_native rows")

# UCSC hg38 wgRna: BED-like table, with miRNA hairpin loci as 0-based half-open intervals.
# GENCODE v49 comprehensive GTF: GRCh38 gene intervals and symbols.
ucsc_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgRna.txt.gz"
gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
ucsc_path <- file.path(ann, "wgRna.txt.gz")
gtf_path <- file.path(ann, "gencode.v49.annotation.gtf.gz")
download_if_missing <- function(url, path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    download.file(url, path, mode = "wb", quiet = FALSE)
  }
  if (file.info(path)$size < 1000) stop("Annotation download looks incomplete: ", path)
}
download_if_missing(ucsc_url, ucsc_path)
download_if_missing(gtf_url, gtf_path)

wg <- read.delim(gzfile(ucsc_path), header = FALSE, sep = "\t", quote = "",
                 comment.char = "", stringsAsFactors = FALSE)
stopifnot(ncol(wg) >= 10)
names(wg)[1:10] <- c("bin", "chrom", "start0", "end", "name", "score", "strand",
                     "thickStart", "thickEnd", "type")
wg <- wg[tolower(trimws(wg$type)) == "mirna", , drop = FALSE]
wg$key <- tolower(wg$name)
wg <- wg[!duplicated(wg$key), , drop = FALSE]

gtf <- read.delim(gzfile(gtf_path), header = FALSE, sep = "\t", quote = "",
                  comment.char = "#", stringsAsFactors = FALSE)
stopifnot(ncol(gtf) == 9)
names(gtf) <- c("chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attr")
gtf <- gtf[gtf$feature == "gene", , drop = FALSE]
attr1 <- function(s, key) {
  m <- regexec(paste0("(?:^|;\\s*)", key, "\\s+\"([^\"]+)\""), s, perl = TRUE)
  v <- regmatches(s, m)
  vapply(v, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
}
genes <- data.frame(chrom = gtf$chrom, start = as.integer(gtf$start), end = as.integer(gtf$end),
                    strand = gtf$strand, gene_id = attr1(gtf$attr, "gene_id"),
                    gene_name = attr1(gtf$attr, "gene_name"), gene_type = attr1(gtf$attr, "gene_type"),
                    stringsAsFactors = FALSE)

mi <- data.frame(miRNA = can$miRNA_names, key = tolower(can$miRNA_names), stringsAsFactors = FALSE)
mi$annotation_status <- ifelse(mi$key %in% wg$key, "mapped_to_UCSC_hg38_wgRna_miRNA_locus", "not_mapped")
mi$chrom <- mi$start0 <- mi$end <- mi$strand <- mi$ucsc_name <- NA
ix <- match(mi$key, wg$key)
ok <- !is.na(ix)
mi$chrom[ok] <- wg$chrom[ix[ok]]
mi$start0[ok] <- wg$start0[ix[ok]]
mi$end[ok] <- wg$end[ix[ok]]
mi$strand[ok] <- wg$strand[ix[ok]]
mi$ucsc_name[ok] <- wg$name[ix[ok]]
mi$identifier_interpretation <- ifelse(ok,
  "input identifier matches a UCSC miRNA locus name; treated as precursor/hairpin locus, not mature 5p/3p",
  "unmapped identifier; precursor identity not resolved from this annotation")

# First retain all fully containing same-strand gene annotations as locus identity context;
# only protein-coding and lncRNA genes are eligible host candidates.
overlap_rows <- list()
for (i in which(ok)) {
  s1 <- as.integer(mi$start0[i]) + 1L
  e1 <- as.integer(mi$end[i])
  hit <- genes$chrom == mi$chrom[i] & genes$strand == mi$strand[i] &
    genes$start <= s1 & genes$end >= e1
  g <- genes[hit, , drop = FALSE]
  if (nrow(g)) {
    g$miRNA <- mi$miRNA[i]
    g$overlap_rule <- "same chromosome and strand; GENCODE gene interval fully contains UCSC precursor interval"
    g$is_host_candidate_type <- !is.na(g$gene_type) & g$gene_type %in% c("protein_coding", "lncRNA")
    g$in_model_3000 <- g$gene_name %in% can$gene_names
    g$available_in_full_mrna <- g$gene_name %in% raw_gene_names
    overlap_rows[[length(overlap_rows) + 1L]] <- g
  }
}
overlapping_genes <- if (length(overlap_rows)) do.call(rbind, overlap_rows) else data.frame()
hosts <- if (nrow(overlapping_genes)) overlapping_genes[overlapping_genes$is_host_candidate_type, , drop = FALSE] else data.frame()
if (nrow(hosts)) hosts$mapping_rule <- "same chromosome and strand; fully contained by GENCODE v49 protein_coding or lncRNA gene interval"
if (nrow(overlapping_genes)) {
  overlapping_genes <- overlapping_genes[order(overlapping_genes$miRNA, overlapping_genes$gene_type,
                                                overlapping_genes$gene_name, overlapping_genes$gene_id), , drop = FALSE]
}
if (nrow(hosts)) {
  hosts <- hosts[order(hosts$miRNA, hosts$gene_name, hosts$gene_id), , drop = FALSE]
  hosts$host_expression_correlation <- NA_real_
  hosts$host_expression_p <- NA_real_
  hosts$host_expression_n_pair <- 0L
}

host_associations <- function(host_df, expression_matrix, analysis_scope) {
  d <- host_df[host_df$gene_name %in% colnames(expression_matrix), , drop = FALSE]
  d$analysis_scope <- analysis_scope
  d$host_expression_correlation <- NA_real_
  d$host_expression_p <- NA_real_
  d$host_expression_n_pair <- 0L
  for (j in seq_len(nrow(d))) {
    x <- expression_matrix[, d$gene_name[j]]
    y <- all_y_native[, d$miRNA[j]]
    pair_ok <- is.finite(x) & is.finite(y)
    n_pair <- sum(pair_ok)
    d$host_expression_n_pair[j] <- n_pair
    if (n_pair >= 3L && sd(x[pair_ok]) > 0 && sd(y[pair_ok]) > 0) {
      d$host_expression_correlation[j] <- suppressWarnings(cor(x[pair_ok], y[pair_ok], method = "spearman"))
      d$host_expression_p[j] <- suppressWarnings(cor.test(x[pair_ok], y[pair_ok], method = "spearman", exact = FALSE)$p.value)
    }
  }
  d <- d[is.finite(d$host_expression_correlation) & is.finite(d$host_expression_p), , drop = FALSE]
  d$p_BH <- p.adjust(d$host_expression_p, method = "BH")
  d
}
association_all_mrna <- host_associations(hosts, all_mrna, "all_18204_mRNA")
association_model_inputs <- host_associations(hosts[hosts$in_model_3000, , drop = FALSE],
  model_input_values, "model_inputs_only_3000")

score <- scores[match(mi$miRNA, scores$miRNA), , drop = FALSE]
expr <- rbind(can$Y_train_native, can$Y_cal_native, can$Y_test_native)
ids <- rownames(expr)
expr <- expr[ids, can$miRNA_names, drop = FALSE]
mean_native <- colMeans(expr, na.rm = TRUE)
var_native <- apply(expr, 2, var, na.rm = TRUE)
summary <- data.frame(miRNA = can$miRNA_names,
  identity_class = ifelse(can$miRNA_names %in% wg$name,
                          "precursor/hairpin locus (UCSC wgRna miRNA)", "precursor-style input name; UCSC locus unmapped"),
  native_expression_mean_366 = mean_native[can$miRNA_names],
  native_expression_variance_366 = var_native[can$miRNA_names],
  expression_scale_note = "miRNA input log2(RPM+1); supported by GDC case comparisons; training z-score is a separate transformation",
  stringsAsFactors = FALSE)
summary <- merge(summary, score[, c("miRNA", "R2", "Pearson", "RMSE", "Rel")], by = "miRNA", sort = FALSE)
summary <- merge(summary, mi[, c("miRNA", "annotation_status", "chrom", "start0", "end", "strand", "ucsc_name")], by = "miRNA", sort = FALSE)
summary$host_gene_candidates <- vapply(summary$miRNA, function(m) {
  if (!nrow(hosts)) return("")
  paste(unique(hosts$gene_name[hosts$miRNA == m]), collapse = ";")
}, character(1))
summary$n_host_gene_candidates <- vapply(summary$miRNA, function(m) {
  if (!nrow(hosts)) return(0L)
  length(unique(hosts$gene_id[hosts$miRNA == m]))
}, integer(1))
summary$host_expression_available_full_mRNA <- vapply(summary$miRNA, function(m) {
  if (!nrow(hosts)) return(FALSE)
  any(hosts$miRNA == m & hosts$available_in_full_mrna)
}, logical(1))
summary$host_gene_candidates_model_inputs_only <- vapply(summary$miRNA, function(m) {
  if (!nrow(hosts)) return("")
  paste(unique(hosts$gene_name[hosts$miRNA == m & hosts$in_model_3000]), collapse = ";")
}, character(1))
summary$n_host_gene_candidates_in_model_3000 <- vapply(summary$miRNA, function(m) {
  if (!nrow(hosts)) return(0L)
  length(unique(hosts$gene_id[hosts$miRNA == m & hosts$in_model_3000]))
}, integer(1))
summary$host_expression_available_model_inputs_only <- summary$n_host_gene_candidates_in_model_3000 > 0L
write.csv(summary, file.path(out, "E5_miRNA_identity_expression_score.csv"), row.names = FALSE, na = "")
write.csv(mi, file.path(out, "E5_UCSC_hairpin_locus_mapping_all846.csv"), row.names = FALSE, na = "")
write.csv(hosts, file.path(out, "E5_host_gene_candidates.csv"), row.names = FALSE, na = "")
write.csv(overlapping_genes, file.path(out, "E5_locus_gene_identity_all_overlaps.csv"), row.names = FALSE, na = "")

# Example selection uses calibration-derived scores only: highest, closest to 0.50, lowest.
ord <- order(score$Rel, decreasing = TRUE)
example_ids <- c(score$miRNA[ord[1]], score$miRNA[which.min(abs(score$Rel - .5))], score$miRNA[ord[length(ord)]])
example_role <- c("high_score", "mid_score_near_0.50", "low_score")
names(example_role) <- example_ids
pred_native <- rbind(
  sweep(sweep(can$pred_cal_z, 2, can$y_scale, "*"), 2, can$y_center, "+"),
  sweep(sweep(can$pred_test_z, 2, can$y_scale, "*"), 2, can$y_center, "+"))
truth_native <- rbind(can$Y_cal_native, can$Y_test_native)
pred_native <- pred_native[rownames(truth_native), , drop = FALSE]
example_rows <- do.call(rbind, lapply(example_ids, function(m) data.frame(
  patient_id = rownames(truth_native), miRNA = m, example_group = unname(example_role[m]),
  Rel = score$Rel[match(m, score$miRNA)], split = c(rep("calibration", nrow(can$Y_cal_native)), rep("test", nrow(can$Y_test_native))),
  observed_native = truth_native[, m], predicted_native = pred_native[, m],
  scale_note = "miRNA input log2(RPM+1); supported by GDC case comparisons; training z-score is a separate transformation")))
write.csv(example_rows, file.path(out, "E5_high_mid_low_example_scatter_data.csv"), row.names = FALSE)
write.csv(data.frame(role = example_role, miRNA = example_ids,
                     Rel = score$Rel[match(example_ids, score$miRNA)]),
          file.path(out, "E5_representative_examples.csv"), row.names = FALSE)

# Describe score/abundance relationships across miRNA features, with uncertainty resampled
# at the patient level to retain cross-miRNA dependence.
score_expr <- merge(scores[, c("miRNA", "Rel")], summary[, c("miRNA", "native_expression_mean_366", "native_expression_variance_366")], by = "miRNA")
relation_point <- c(
  mean = cor(score_expr$Rel, score_expr$native_expression_mean_366, method = "spearman"),
  variance = cor(score_expr$Rel, score_expr$native_expression_variance_366, method = "spearman"))
bootstrap_B <- as.integer(Sys.getenv("E5_BOOTSTRAP_REPLICATES", "2000"))
if (is.na(bootstrap_B) || bootstrap_B < 1L) stop("E5_BOOTSTRAP_REPLICATES must be a positive integer")
set.seed(20261014)
boot_rho <- matrix(NA_real_, nrow = bootstrap_B, ncol = 2L,
                   dimnames = list(NULL, c("mean", "variance")))
clamp01 <- function(z) pmax(0, pmin(1, z))
for (b in seq_len(bootstrap_B)) {
  itrain <- sample.int(nrow(can$Y_train_native), nrow(can$Y_train_native), replace = TRUE)
  ical <- sample.int(nrow(can$Y_cal_native), nrow(can$Y_cal_native), replace = TRUE)
  itest <- sample.int(nrow(can$Y_test_native), nrow(can$Y_test_native), replace = TRUE)
  yboot <- rbind(can$Y_train_native[itrain, , drop = FALSE],
                 can$Y_cal_native[ical, , drop = FALSE],
                 can$Y_test_native[itest, , drop = FALSE])
  mu <- colMeans(yboot, na.rm = TRUE)
  variance <- colSums(sweep(yboot, 2, mu, "-")^2, na.rm = TRUE) / (nrow(yboot) - 1)
  yt <- can$Y_cal_z[ical, , drop = FALSE]
  yp <- can$pred_cal_z[ical, , drop = FALSE]
  r2 <- 1 - colSums((yt - yp)^2, na.rm = TRUE) /
    colSums(sweep(yt, 2, colMeans(yt, na.rm = TRUE), "-")^2, na.rm = TRUE)
  yc <- sweep(yt, 2, colMeans(yt, na.rm = TRUE), "-")
  pc <- sweep(yp, 2, colMeans(yp, na.rm = TRUE), "-")
  pearson <- colSums(yc * pc, na.rm = TRUE) /
    sqrt(colSums(yc^2, na.rm = TRUE) * colSums(pc^2, na.rm = TRUE))
  rel <- .5 * clamp01(r2) + .5 * clamp01(pearson)
  rel[!is.finite(rel)] <- 0
  boot_rho[b, "mean"] <- suppressWarnings(cor(rel, mu, method = "spearman", use = "complete.obs"))
  boot_rho[b, "variance"] <- suppressWarnings(cor(rel, variance, method = "spearman", use = "complete.obs"))
}
bootstrap_result <- do.call(rbind, lapply(c("mean", "variance"), function(k) {
  z <- boot_rho[, k]
  z <- z[is.finite(z)]
  data.frame(expression_statistic = if (k == "mean") "native_expression_mean" else "native_expression_variance",
    point_rho = unname(relation_point[k]), bootstrap_median = if (length(z)) median(z) else NA_real_,
    ci_2.5 = if (length(z)) unname(quantile(z, .025, names = FALSE)) else NA_real_,
    ci_97.5 = if (length(z)) unname(quantile(z, .975, names = FALSE)) else NA_real_,
    B_valid = length(z), B_requested = bootstrap_B,
    resampling_unit = "patients sampled with replacement separately within train (222), calibration (72), and test (72); native expression summaries recomputed over the combined 366 draws; calibration truth and fixed-model predictions resampled together to recompute Rel" )
}))
write.csv(bootstrap_result, file.path(out, "E5_score_expression_patient_bootstrap.csv"), row.names = FALSE)

# Save patient-level aligned host associations for both the full matrix and model inputs.
write.csv(association_all_mrna, file.path(out, "E5_host_expression_associations_all_mRNA.csv"), row.names = FALSE, na = "")
write.csv(association_model_inputs, file.path(out, "E5_host_expression_associations_model_inputs_only.csv"), row.names = FALSE, na = "")

matched_label_test <- function(association, scope, seed = 20261014L, B = 10000L) {
  matched <- summary[is.finite(summary$Rel) & is.finite(summary$native_expression_mean_366) &
                       is.finite(summary$native_expression_variance_366), , drop = FALSE]
  matched$score_group <- ifelse(matched$Rel >= .7, "high", "background")
  matched$mean_bin <- pmin(10L, ceiling(rank(matched$native_expression_mean_366, ties.method = "average") / nrow(matched) * 10))
  matched$variance_bin <- pmin(10L, ceiling(rank(log1p(matched$native_expression_variance_366), ties.method = "average") / nrow(matched) * 10))
  matched$stratum <- paste(matched$mean_bin, matched$variance_bin, sep = "_")
  assoc_abs <- aggregate(abs(host_expression_correlation) ~ miRNA, association, mean)
  names(assoc_abs)[2] <- "mean_abs_host_rho"
  matched <- merge(matched, assoc_abs, by = "miRNA", all.x = TRUE)
  high_all <- matched[matched$score_group == "high" & is.finite(matched$mean_abs_host_rho), , drop = FALSE]
  background_all <- matched[matched$score_group == "background" & is.finite(matched$mean_abs_host_rho), , drop = FALSE]
  common_strata <- intersect(unique(high_all$stratum), unique(background_all$stratum))
  common <- matched[matched$stratum %in% common_strata & is.finite(matched$mean_abs_host_rho), , drop = FALSE]
  strata <- sort(unique(common$stratum))
  nh_common <- sum(common$score_group == "high")
  nb_common <- sum(common$score_group == "background")
  observed <- p <- NA_real_
  B_run <- 0L
  if (length(strata) > 0L && nh_common > 0L && nb_common > 0L) {
    weights <- vapply(strata, function(st) sum(common$stratum == st & common$score_group == "high") / nh_common, numeric(1))
    observed_by_stratum <- vapply(strata, function(st) {
      z <- common[common$stratum == st, , drop = FALSE]
      mean(z$mean_abs_host_rho[z$score_group == "high"]) - mean(z$mean_abs_host_rho[z$score_group == "background"])
    }, numeric(1))
    observed <- sum(weights * observed_by_stratum)
    set.seed(seed)
    null <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      perm_by_stratum <- vapply(strata, function(st) {
        z <- common[common$stratum == st, , drop = FALSE]
        nh <- sum(z$score_group == "high")
        hi <- sample.int(nrow(z), size = nh, replace = FALSE)
        mean(z$mean_abs_host_rho[hi]) - mean(z$mean_abs_host_rho[-hi])
      }, numeric(1))
      null[b] <- sum(weights * perm_by_stratum)
    }
    p <- (1 + sum(abs(null) >= abs(observed))) / (B + 1)
    B_run <- B
  }
  data.frame(analysis_scope = scope, n_high_with_host_association = nrow(high_all),
    n_background_with_host_association = nrow(background_all), n_common_support_strata = length(strata),
    n_high_common_support = nh_common, n_background_common_support = nb_common,
    n_high_without_common_support = nrow(high_all) - nh_common,
    observed_stratum_weighted_mean_abs_rho_difference = observed,
    label_permutation_p_two_sided = p, permutations = B_run,
    exchangeability_rule = "within each joint expression-mean and log1p-variance decile, randomly relabel high/background without replacement while preserving the observed high count; only strata containing both labels enter; effect is high-count-weighted mean of within-stratum high minus background mean absolute Spearman rho",
    interpretation = "exploratory co-expression analysis; exchangeability is conditional on these measured strata; sparse common support limits precision; not causal evidence")
}
matched_all <- matched_label_test(association_all_mrna, "all_18204_mRNA")
matched_model <- matched_label_test(association_model_inputs, "model_inputs_only_3000")
write.csv(rbind(matched_all, matched_model), file.path(out, "E5_host_expression_stratified_label_permutation.csv"), row.names = FALSE)

meta <- list(annotation_sources = list(
  list(name = "UCSC hg38 wgRna table", url = ucsc_url, assembly = "GRCh38/hg38", interpretation = "miRNA locus intervals and names"),
  list(name = "GENCODE human v49 comprehensive annotation", url = gtf_url, assembly = "GRCh38", interpretation = "gene intervals, gene symbols, and gene types")),
  annotation_file_md5 = as.list(tools::md5sum(c(ucsc_path, gtf_path))),
  mapping_rule = "case-insensitive exact input identifier to UCSC miRNA locus name; precursor interval must be fully contained in same-strand GENCODE gene interval for host candidate",
  access_date = format(Sys.Date(), "%Y-%m-%d"), n_input = length(can$miRNA_names), n_locus_mapped = sum(ok),
  host_gene_type_rule = "host candidates restricted to GENCODE gene_type protein_coding or lncRNA; overlapping miRNA gene records and other types are retained only in E5_locus_gene_identity_all_overlaps.csv, not counted as host genes",
  n_with_host_candidate = if (nrow(hosts)) length(unique(hosts$miRNA)) else 0L,
  n_host_candidate_rows = nrow(hosts),
  host_candidate_gene_type_counts = if (nrow(hosts)) as.list(setNames(as.integer(table(hosts$gene_type)), names(table(hosts$gene_type)))) else list(),
  n_overlapping_gene_records_all_types = nrow(overlapping_genes),
  n_overlapping_miRNA_gene_identity_records_excluded_from_hosts = if (nrow(overlapping_genes)) sum(overlapping_genes$gene_type == "miRNA", na.rm = TRUE) else 0L,
  n_with_host_candidate_among_Rel_ge_0.7 = if (nrow(hosts)) length(intersect(unique(hosts$miRNA), scores$miRNA[scores$Rel >= 0.7])) else 0L,
  n_with_host_expression_feature_all_18204_mRNA = if (nrow(hosts)) length(unique(hosts$miRNA[hosts$available_in_full_mrna])) else 0L,
  n_with_host_expression_feature_among_Rel_ge_0.7_all_18204_mRNA = if (nrow(hosts)) length(intersect(unique(hosts$miRNA[hosts$available_in_full_mrna]), scores$miRNA[scores$Rel >= 0.7])) else 0L,
  n_with_host_expression_feature_model_inputs_only_3000 = if (nrow(hosts)) length(unique(hosts$miRNA[hosts$in_model_3000])) else 0L,
  n_with_host_expression_feature_among_Rel_ge_0.7_model_inputs_only_3000 = if (nrow(hosts)) length(intersect(unique(hosts$miRNA[hosts$in_model_3000]), scores$miRNA[scores$Rel >= 0.7])) else 0L,
  raw_mrna_patient_alignment = raw_mrna_alignment,
  n_host_expression_pairs_all_18204_mRNA = nrow(association_all_mrna),
  n_host_expression_pairs_model_inputs_only_3000 = nrow(association_model_inputs),
  n_host_expression_pairs_BH_q_below_0_05_all_18204_mRNA = sum(association_all_mrna$p_BH < .05, na.rm = TRUE),
  n_host_expression_pairs_BH_q_below_0_05_model_inputs_only_3000 = sum(association_model_inputs$p_BH < .05, na.rm = TRUE),
  score_expression_bootstrap_replicates_requested = bootstrap_B,
  score_expression_bootstrap_replicates_valid_mean = bootstrap_result$B_valid[bootstrap_result$expression_statistic == "native_expression_mean"],
  score_expression_bootstrap_replicates_valid_variance = bootstrap_result$B_valid[bootstrap_result$expression_statistic == "native_expression_variance"],
  score_expression_bootstrap_seed = 20261014,
  expression_unit = "miRNA log2(RPM+1) supported by five GDC cases; mRNA log2(RSEM norm_count+1) verified in all primary cells",
  limitations = c("UCSC locus identity does not identify mature 5p/3p products", "no validated mature-target database was queried or inferred", "same-strand containment is a candidate host relation, not transcriptional proof", "host-gene expression correlations use the same 366 patients and are observational", "no independent tumor-purity or cell-composition field available in canonical input"))
jsonlite::write_json(meta, file.path(out, "E5_annotation_metadata.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
cat("E5 complete: mapped", sum(ok), "/", length(ok), "precursor-style identifiers; host candidate loci",
    meta$n_with_host_candidate, "; full-mRNA host loci", meta$n_with_host_expression_feature_all_18204_mRNA,
    "; model-input host loci", meta$n_with_host_expression_feature_model_inputs_only_3000,
    "; bootstrap B", bootstrap_B, "\n")
