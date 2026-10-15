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
| `verification/` | Scripts that recompute the reported revision statistics directly from result tables, used to check every number quoted in the response letter |
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
`public_compare_xena_hiseqv2_366.R` (366 x 18,204 cells identical) and the verified ASA units are
documented in the manuscript.

## How to run

The scripts use repository-root-relative paths anchored at `code/`, `new/code/` and
`new/results/`, i.e. the layout below. To execute them unchanged, recreate this layout and run
everything from the repository root:

```
repo/
  data/data.rda
  code/            <- files from OA/main_pipeline
  new/code/        <- files from OA/revision_analysis, external_data, figures, verification
  new/results/     <- outputs (created automatically)
```

Main pipeline:

```
Rscript code/StrongMiRRecon_v2.R
Rscript code/run_external_generalization.R --project-dir . --external-data <cohort>.rda
```

Revision pipeline (order matters; long-running steps can be checkpointed and restarted safely):

```
Rscript new/code/prepare_canonical.R                  # E0 canonical fixed-split input
Rscript new/code/data_description.R                   # E0 descriptive tables
Rscript new/code/repeated_training.R all              # E1  20 patient-level splits
Rscript new/code/repeated_evaluation.R watch          # E1  per-split evaluation, auto-resumes
Rscript new/code/shuffled_baseline.R all              # 20 training-pairing permutations
Rscript new/code/downstream_analysis.R <canonical.rds> new/results/E3   # E3 fixed split
Rscript new/code/score_external_analysis.R --canonical <canonical.rds> --section E2 --bootstrap 2000 --output-root new/results
Rscript new/code/uncertainty_completion.R --bootstrap 2000 --section models
Rscript new/code/uncertainty_completion.R --bootstrap 2000 --section scores
Rscript new/code/uncertainty_completion.R --bootstrap 2000 --section sample
Rscript new/code/build_chol_harmonized_input.R        # external input, same RSEM pipeline
Rscript new/code/score_external_analysis.R --section E4 --bootstrap 2000 --external-input <chol.rda> --output-root new/results/harmonized
Rscript new/code/build_multicohort_inputs.R           # 7 additional TCGA cohorts
Rscript new/code/run_multicohort_e4.ps1               # per-cohort transfer evaluation
Rscript new/code/summarize_multicohort_e4.R           # 8-cohort summary
Rscript new/code/incremental_prognosis_comparison.R   # three-model prognostic comparison
Rscript new/code/biological_annotation.R              # E5 identity and host analysis
Rscript new/code/gen_fig_*.R                          # figure PDF/PNG
```

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
