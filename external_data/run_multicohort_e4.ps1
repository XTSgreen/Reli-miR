# Run the frozen E4 transfer analysis for each TCGA cohort.
# Requires inputs from build_multicohort_inputs.R.
$ErrorActionPreference = 'Continue'
$cohorts = 'KIRC','LUAD','COAD','STAD','PRAD','THCA','HNSC'
foreach ($co in $cohorts) {
  $rda = "new\results\E4_multicohort\inputs\$co\external_${co}_matched_RSEM.rda"
  if (-not (Test-Path $rda)) { Write-Host "$co input missing - skipped"; continue }
  $done = "new\results\harmonized_multicohort\$co\E4\E4_external_overall_metrics.csv"
  if (Test-Path $done) { Write-Host "$co already analyzed - skipped"; continue }
  Write-Host "=== E4 transfer: $co ==="
  & Rscript 'new/code/score_external_analysis.R' --section E4 --bootstrap 2000 `
    --external-input $rda --output-root "new/results/harmonized_multicohort/$co"
  Write-Host "=== exit: $co exitcode=$LASTEXITCODE ==="
}
Write-Host 'ALL E4 TRANSFER RUNS DONE'
