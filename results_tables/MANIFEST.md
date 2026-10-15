# Result table manifest

Every number quoted in the manuscript response letter and in the revision can be checked against
the tables in this folder. Run `Rscript verification/verify_all.R` from the repository root to
recompute all of them; the script exits with status 1 if any published number does not reproduce.

Raw intermediate outputs are **not** redistributed here. Two multi-megabyte per-miRNA frequency
tables (`E1_miRNA_selection_frequency.csv`, ~5 MB; `E3_winner_by_seed.csv`, ~4.7 MB) are excluded;
they are available from the authors on request.

## E0 — data provenance, units and training settings

| Table | Content | Produced by | Supports |
|---|---|---|---|
| `E0_audit.json` | Version and provenance audit of the frozen analysis | `revision_analysis/prepare_canonical.R` | E0 version reconciliation |
| `E0_model_input_gene_list.csv` | The 3,000 model input genes | `revision_analysis/prepare_canonical.R` | Reviewer 3 minor point on mRNA features |
| `E0_miRNA_input_descriptive_statistics.csv` | Per-miRNA input descriptives for all 846 miRNAs | `revision_analysis/data_description.R` | Per-transcript statistics requested by Reviewer 3 |
| `E0_fixed_split_patient_manifest.csv` | Patient assignment of the fixed 222/72/72 split | `revision_analysis/prepare_canonical.R` | Reproducibility of the submitted split |
| `E0_clinical_group_counts.csv` | Stage and event counts | `revision_analysis/data_description.R` | Downstream group sizes |
| `E0_original_MLP_training_settings.csv` | Recovered settings of the original R MLP | `revision_analysis/prepare_canonical.R` | Neural-network detail requested by Reviewer 4 |
| `public_provenance_search.md` | Public-source verification narrative | `external_data/public_provenance_verify.R` | Reviewer 3/4 data-source points |
| `public_LIHC_Xena_HiSeqV2_match_summary_366.csv` | Cell-level match of the 366 x 18,204 primary matrix against Xena HiSeqV2 | `external_data/public_compare_xena_hiseqv2_366.R` | mRNA platform and normalisation |
| `public_primary_OS_match.csv` | Overall-survival time match against GDC | `external_data/public_OS_compare.R` | Time unit (days) |
| `public_gdc_sample_value_match.csv` | miRNA value match for sampled GDC cases | `external_data/public_gdc_sample_compare.R` | miRNA unit log2(RPM+1) |

## E1 — 20-split repeated training and evaluation

| Table | Content | Produced by | Supports |
|---|---|---|---|
| `E1_evaluation_progress.csv` | Completion state of the 20 evaluations | `revision_analysis/repeated_evaluation.R` | Protocol completion |
| `E1_final_model_performance_distribution.csv` | Mean/median/quantiles of test metrics over 20 splits, per model | `revision_analysis/repeated_evaluation.R` | Model comparison over splits |
| `E1_winner_frequency_partial.csv` | Training-internal selected-model counts (14 / 6 / 0) | `revision_analysis/repeated_evaluation.R` | Selection instability |
| `E1_model_performance_by_seed.csv` | Per-split metrics, ranks and selected hyper-parameters | `revision_analysis/repeated_evaluation.R` | Paired differences |
| `E1_selected_hyperparameters_by_seed.csv` | Chosen k, alpha, lambda index and ratio per split | `revision_analysis/repeated_evaluation.R` | Hyper-parameter reproducibility |
| `E1_hyperparameter_selection_frequency.csv` | Frequency of each configuration | `revision_analysis/repeated_evaluation.R` | Configuration stability |
| `E1_original_33_selection_frequency.csv` | Per-miRNA retention frequency of the original 33 | `revision_analysis/repeated_evaluation.R` | Reviewer 4 stability question |
| `E1_Rel_threshold_selection_frequency.csv` | Retention frequency per miRNA, model and threshold | `revision_analysis/repeated_evaluation.R` | Threshold lists |
| `E1_selection_quality_by_seed.csv` | Selected-subset test quality per method and coverage | `revision_analysis/repeated_evaluation.R` | Equal-coverage score comparison |
| `E1_threshold_counts_by_seed.csv` | Retained counts at Rel >= 0.5/0.6/0.7 per split and model | derived from per-split evaluations | Threshold count distribution |
| `E1_score_associations_by_seed.csv` | Quality-score vs test-error Spearman per method and split | derived from per-split evaluations | Score ablation over splits |
| `E1_toplist_jaccard_summary.csv` | Equal-size top-list Jaccard summary | derived from per-split evaluations | List stability |
| `E1_Rel_threshold_pairwise_Jaccard.csv` | Pairwise Jaccard of Rel-threshold lists | `revision_analysis/repeated_evaluation.R` | List stability |
| `E1_shuffled/E1_shuffled_20_permutation_summary.csv` | 20 training-pairing permutations | `revision_analysis/shuffled_baseline.R` | Negative control (Reviewer 4) |
| `E1_shuffled/README_shuffled.md` | Notes on the realised Ridge baseline | `revision_analysis/shuffled_baseline.R` | Baseline naming correction |

## E2 — score ablation, thresholds, sample level

| Table | Content | Produced by | Supports |
|---|---|---|---|
| `E2_complete_model_metrics_point.csv` | Fixed-split point estimates | `revision_analysis/score_external_analysis.R` | Model comparison |
| `E2_complete_model_bootstrap_summary.csv` | Paired patient bootstrap intervals | `revision_analysis/uncertainty_completion.R` | Reviewer 1/4 uncertainty request |
| `E2_complete_model_paired_contrasts.csv` | Model differences with intervals | `revision_analysis/uncertainty_completion.R` | Small-difference claim |
| `E2_complete_score_point_reproduction_check.csv` | Point Spearman of ten scoring rules vs test RMSE | `revision_analysis/score_external_analysis.R` | Score ablation |
| `E2_complete_score_joint_bootstrap_summary.csv` | Joint train/calibration/test resampling intervals | `revision_analysis/uncertainty_completion.R` | Score ablation intervals |
| `E2_complete_sample_bootstrap_summary.csv` | Sample-level score components with intervals | `revision_analysis/uncertainty_completion.R` | Reviewer 1 sample-level point |
| `E2_Rel_threshold_coverage_FNR_point.csv` | Precision, recall, FNR and coverage over thresholds 0.10-0.90 | `revision_analysis/score_external_analysis.R` | Reviewer 3 false-negative question |
| `E2_FNR_bootstrap_summary.csv` | Bootstrap intervals for FNR | `revision_analysis/uncertainty_completion.R` | Threshold cost |
| `E2_FNR_label_sensitivity_point.csv` | Alternative operational-label definitions | `revision_analysis/score_external_analysis.R` | Label sensitivity |
| `E2_equal_coverage_score_comparison.csv` | Same-coverage comparison of scoring rules | `revision_analysis/score_external_analysis.R` | Reviewer 2 equal-coverage point |
| `README_E2.md` | Section notes | — | — |

## E3 — downstream preservation and size-matched controls

| Table | Content | Produced by | Supports |
|---|---|---|---|
| `E3_subset_summary_named.csv` | Fixed-split downstream metrics for named subsets (random rows removed) | `revision_analysis/downstream_analysis.R` | Table 3 of the manuscript |
| `E3_core_patient_bootstrap_CI_fast_B1000.csv` | Patient bootstrap intervals (B=1000) for core metrics | `revision_analysis/downstream_analysis.R` | Reviewer 4 interval request |
| `E3_paired_comparison_bootstrap_CI.csv` | Rel versus random and versus simple-score contrasts | `revision_analysis/downstream_paired_comparison.R` | Size-matched control claim |
| `E3_cluster_ARI_NMI_named.csv` | ARI/NMI for k=2-6, named subsets | `revision_analysis/downstream_analysis.R` | Reviewer 4 Table S4 point |
| `E3_downstream_by_split.csv` | Per-split values for Rel, Pearson, calibration-RMSE and OOF-RMSE subsets | derived from per-split evaluations | Control consistency over splits |
| `E3_random_mean_by_split.csv` | Per-split random-list means | derived from per-split evaluations | Control consistency over splits |
| `E3_random_subset_distribution_summary.csv` | Fixed-split random-subset distribution summary (1,000 lists per size) | derived from `E3_random_subset_distribution.csv` | Same-size control requirement |
| `README_E3.md` | Section notes | — | — |

## E4 — external transfer (CHOL and eight cohorts)

| Table | Content | Produced by | Supports |
|---|---|---|---|
| `E4_CHOL/E4_external_overall_metrics.csv` | Overall external metrics, 742 shared miRNAs | `revision_analysis/score_external_analysis.R` | External evaluation |
| `E4_CHOL/E4_complete_external_bootstrap_summary.csv` | Patient bootstrap (B=2000) point estimates and intervals | `revision_analysis/uncertainty_completion.R` | Reviewer 1/2 external intervals |
| `E4_CHOL/E4_external_patient_bootstrap_summary.csv` | Bootstrap summary in the section-local format | `revision_analysis/score_external_analysis.R` | Cross-check of the above |
| `E4_CHOL/E4_external_metrics_by_miRNA.csv` | Per-miRNA external metrics and frozen internal score | `revision_analysis/score_external_analysis.R` | Reviewer 4 external reliability transfer |
| `E4_CHOL/E4_Rel_threshold_external_coverage.csv` | External coverage and performance of frozen internal lists | `revision_analysis/score_external_analysis.R` | High-reliability subset in the external cohort |
| `E4_CHOL/E4_equal_size_selection_summary.csv` | Equal-size external selection comparison | `revision_analysis/score_external_analysis.R` | Size-matched external control |
| `E4_CHOL/E4_run_metadata.csv` | Run metadata and coverage | `revision_analysis/score_external_analysis.R` | Provenance |
| `E4_CHOL/README_E4_harmonized.md` | Section notes | `external_data/finalize_harmonized_e4.py` | — |
| `E4_multicohort/multicohort_E4_summary.csv` | Eight-cohort transfer summary (3,269 patients) | `revision_analysis/summarize_multicohort_e4.R` | Cross-cohort ranking transfer |
| `E4_multicohort/README_E4_multicohort.md` | Multi-cohort section notes | `external_data/finalize_harmonized_e4.py` | — |
| `E4_multicohort/PREREGISTRATION_multicohort_E4.md` | Frozen design written before any evaluation | — | Defensive record against post-hoc selection |
| `E4_multicohort/inputs/input_build_summary.csv` | Per-cohort input build summary | `external_data/build_multicohort_inputs.R` | Input coverage 3000/3000 and 846/846 |

## E5 — biological annotation and prognostic comparison

| Table | Content | Produced by | Supports |
|---|---|---|---|
| `E5_miRNA_identity_expression_score.csv` | Identity class, location, candidate hosts and expression features for all 846 miRNAs | `revision_analysis/biological_annotation.R` | Editor and Reviewer 3 identity question |
| `E5_host_gene_candidates.csv` | Candidate host relationships (765 rows, 724 miRNAs) | `revision_analysis/biological_annotation.R` | Host-gene analysis |
| `E5_host_expression_stratified_label_permutation.csv` | Stratified label-permutation test of host co-expression | `revision_analysis/biological_annotation.R` | Biological-null result |
| `E5_host_expression_associations_model_inputs_only.csv` | Host associations within the model input range | `revision_analysis/biological_annotation.R` | Non-detection within input range |
| `E5_threemodel_survival_comparison.csv` | Three-model prognostic comparison with B=2000 intervals | `revision_analysis/incremental_prognosis_comparison.R` | Reviewer 4 independent-information point |
| `E5_threemodel_run_metadata.csv` | Selected lambdas, active features and event counts | `revision_analysis/incremental_prognosis_comparison.R` | Provenance |
| `E5/PREREGISTRATION_threemodel_comparison.md` | Frozen design and pre-run amendments | — | Defensive record |
| `E5/README_E5.md` | Section notes | — | — |

## Reading conventions

A reliability threshold of 0.7 is stored with floating-point residue
(`0.7000000000000001`) in some tables; comparisons in the verification script use a tolerance.
Retention-frequency tables are reported per miRNA and per split, so a miRNA absent from a split's
list has frequency zero rather than a missing value.
