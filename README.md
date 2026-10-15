# Reli-miR — Reliability-Aware Virtual miRNA Reconstruction from mRNA

Code accompanying the manuscript *Reli-miR: A Reliability-Aware Framework for Virtual miRNA
Reconstruction from mRNA Expression with Downstream Preservation Assessment* and its revision.

The framework reconstructs miRNA expression from paired mRNA-miRNA transcriptomes and, more
importantly, returns a **per-miRNA reliability score** that is estimated on held-out data and used
to select a trustworthy subset for downstream analyses. The revision adds a 20-split repeated
training/evaluation protocol, a permutation negative control, equal-size random-subset controls,
an eight-cohort external transfer analysis, and a three-model prognostic comparison.

## Repository layout

| Folder | Contents |
|---|---|
| `main_pipeline/` | Submitted pipeline: the full sequential analysis script `StrongMiRRecon_v2.R`, external generalization, PCA_ElasticNet and Pathway_Ridge reruns, GPU (PyTorch) MLP/Transformer comparators |
| `revision_analysis/` | Revision analyses: canonical input build, 20-split repeated training and evaluation, downstream preservation and paired comparison, score ablation and external scoring, uncertainty completion (bootstrap), shuffled negative control, biological annotation, three-model prognostic comparison, multi-cohort aggregation |
| `external_data/` | Public-source provenance verification (GDC/Xena), harmonized CHOL input construction, seven-cohort TCGA input construction, cohort identification, and GDC miRNA bulk-download helper scripts |
| `figures/` | Shared style helper and per-figure generation scripts (E1 stability, E2 score ablation and coverage, E3 equal-size controls, E4 transfer and multi-cohort, E5 expression associations) |
| `results_tables/` | Curated result tables behind every reported number, with per-table provenance in `results_tables/MANIFEST.md`. Raw intermediate outputs and the two multi-megabyte frequency tables are not redistributed |
| `verification/` | Single verification entry point `verify_all.R`, which recomputes every published statistic from `results_tables/` and exits with a non-zero status if any check fails |
| `environment/` | `sessionInfo` and run-environment records (R version and package versions) |

## Requirements

* R 4.6.1 with `glmnet`, `survival`, `pls`, `nnet`, `msigdbr`, `readr`, `tidyverse`
* Python 3 with `python-docx` and `pandas` (document extraction and GPU-input export only)
* PyTorch with CUDA (supplementary GPU comparators only; not required for the main results)

Exact versions are recorded in `environment/33_sessionInfo.txt` and `environment/run_environment.txt`.

## Data

The scripts expect paired mRNA-miRNA data in an R workspace with three objects:

```
clinical : patient_id, OS_time, OS_status, age, gender, stage
mrna     : patient_id + mRNA expression matrix (log2(RSUM/normalised count + 1))
mirna    : patient_id + miRNA expression matrix (log2(RPM + 1), precursor/hairpin identifiers)
```

Place it at `data/data.rda`. The primary cohort is TCGA-LIHC (366 primary tumours); the external
cohorts are TCGA-CHOL (35) and seven additional TCGA solid-tumour cohorts. **Raw data are not
redistributed here.** They are publicly available from the GDC Data Portal and the UCSC Xena
public hub; `external_data/` contains the retrieval and harmonisation scripts. One-step
cross-comparison of the local matrices against the official Xena HiSeqV2 release is provided in
`public_compare_xena_hiseqv2_366.R` (366 x 18,204 cells identical) and the verified expression units
are documented in the manuscript. The curated tables behind every reported number are published
under `results_tables/`, with provenance in `results_tables/MANIFEST.md`; raw result intermediates
are not redistributed.

## How to run

The scripts use repository-root-relative paths anchored at `code/`, `new/code/` and
`new/results/`, i.e. the layout below. To execute them unchanged, recreate this layout and run
everything from the repository root:

```
repo/
  data/data.rda
  code/            <- files from OA/main_pipeline
  new/code/        <- files from OA/revision_analysis, external_data, figures
  new/results/     <- outputs (created automatically)
  results_tables/  <- files from OA/results_tables
  verification/    <- files from OA/verification (reads results_tables/)
```

Main pipeline:

```
Rscript code/StrongMiRRecon_v2.R
Rscript code/export_gpu_model_inputs.R                # model input matrices for the GPU comparators
Rscript code/run_external_generalization.R --project-dir . --external-data <cohort>.rda
```

Revision pipeline. Steps must be run in this order from the repository root; the long-running
training and evaluation steps are checkpointed and can be restarted safely. `code/export_gpu_model_inputs.R`
records the model input matrices that the frozen analysis builds on, so it runs before the E0 step.

```
# Stage 0 - model inputs exported by the submitted pipeline
Rscript code/export_gpu_model_inputs.R

# Stage 1 - E0 foundation
Rscript new/code/prepare_canonical.R                  # canonical fixed-split input
Rscript new/code/data_description.R                   # descriptive tables

# Stage 2 - E1 repeated protocol
Rscript new/code/repeated_training.R all              # 20 patient-level splits
Rscript new/code/repeated_evaluation.R watch          # per-split evaluation, auto-resumes
Rscript new/code/shuffled_baseline.R all              # 20 training-pairing permutations

# Stage 3 - fixed-split E2 and E3
Rscript new/code/downstream_analysis.R <canonical.rds> new/results/E3
Rscript new/code/score_external_analysis.R --canonical <canonical.rds> --section E2 --bootstrap 2000 --output-root new/results
Rscript new/code/uncertainty_completion.R --bootstrap 2000 --section models
Rscript new/code/uncertainty_completion.R --bootstrap 2000 --section scores
Rscript new/code/uncertainty_completion.R --bootstrap 2000 --section sample

# Stage 4 - E4 external cohorts
Rscript new/code/build_chol_harmonized_input.R        # same-RSEM CHOL input
Rscript new/code/verify_chol_harmonized_input.R       # input validation
Rscript new/code/score_external_analysis.R --section E4 --bootstrap 2000 --external-input <chol.rda> --output-root new/results/harmonized
Rscript new/code/uncertainty_completion.R --section external --bootstrap 2000 --external-input <chol.rda> --output-root new/results/harmonized
powershell -ExecutionPolicy Bypass -File new/code/fetch_gdc_mirna_manifests.ps1   # GDC miRNA file manifests
powershell -ExecutionPolicy Bypass -File new/code/download_gdc_mirna.ps1          # bulk miRNA quantification download
Rscript new/code/identify_xena_files.R                # identify each downloaded Xena HiSeqV2 matrix by cohort
Rscript new/code/match_xena_to_project.R              # confirm the cohort assignment against GDC case lists
Rscript new/code/build_multicohort_inputs.R           # seven additional TCGA cohort inputs
powershell -ExecutionPolicy Bypass -File new/code/run_multicohort_e4.ps1          # per-cohort transfer evaluation
Rscript new/code/summarize_multicohort_e4.R           # eight-cohort summary

# Stage 5 - E5
Rscript new/code/incremental_prognosis_comparison.R   # three-model prognostic comparison
Rscript new/code/biological_annotation.R              # identity and host analysis

# Figures
Rscript new/code/gen_fig_*.R                          # figure PDF/PNG

# Verification of the published numbers
Rscript verification/verify_all.R
```

Two notes on Stage 4. The three `.ps1` helpers are PowerShell scripts and must be invoked through
`powershell -ExecutionPolicy Bypass -File`; they are Windows-specific. The Xena HiSeqV2 matrices
themselves are not fetched by these scripts and must be downloaded beforehand from the UCSC Xena
public hub; `identify_xena_files.R` then resolves which cohort each local file belongs to.

## Verifying the published numbers

```
Rscript verification/verify_all.R                 # verifies against results_tables/
Rscript verification/verify_all.R <tables-root>   # alternative location
```

The script recomputes 72 published statistics directly from `results_tables/`, prints the source
file for every check, and reports a PASS/FAIL count. It exits with status 0 when all checks pass
and status 1 when any check fails, so it can be used as a CI gate. `results_tables/MANIFEST.md`
maps each table to the analysis script that produced it and to the manuscript or response-letter
claim it supports.

## Frozen analysis parameters

| Item | Value |
|---|---|
| Repeated splits | 20 patient-level stratified splits, seeds `20261014`–`20261033`; train 222 / calibration 72 / test 72 |
| Model and hyper-parameter selection | moved inside the training set (5-fold CV), calibration set used only for reliability, test set used only for evaluation |
| Lambda protocol | `trainCV_lambda_default_v2`: 100-point path, minimum ratio 0.01 when n < p and 0.0001 when n >= p |
| PCA dimensions / ElasticNet alpha | k in {20, 50, 100, 150, 200}; alpha in {0.1, 0.3, 0.5, 0.7, 0.9} |
| Reliability thresholds | 0.5, 0.6, 0.7 (and a 0.10–0.90 curve at step 0.05) |
| Equal-size coverage | 3.9%, 8.7%, 20.9% of evaluable miRNAs |
| Bootstrap | B = 2000 (external and paired-model intervals; seed `20261017`); B = 1000 (E3 random-subset and patient intervals; seed `20261018`); three-model comparison B = 2000, seed `20261019` |
| Random controls | 1000 equal-size random lists per size on the fixed split; 200 per size in each repeated split (seed `20261018`) |
| Shuffled control | 20 training-pairing permutations; realised baseline is multi-response **Ridge**, not ElasticNet |

## Interpretation boundary

The reported results are deliberately asymmetric and the code is published to make that
checkable: reliability filtering is consistently better than equal-size random subsets, but the
composite score is not consistently better than simpler scores; ranking transfer is directionally
consistent across all eight cohorts while absolute reconstruction is negative in every cohort; and
no incremental prognostic value was established. Reported numbers are recomputed from the result
tables by the scripts in `verification/`.

## License and citation

Released under the MIT License; see `LICENSE`. Citation details follow the manuscript.
