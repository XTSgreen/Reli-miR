# Fetch GDC miRNA quantification file manifests for 7 TCGA cohorts.
$ErrorActionPreference = 'Stop'
$cohorts = @{ 'KIRC' = 'TCGA-KIRC'; 'LUAD' = 'TCGA-LUAD'; 'COAD' = 'TCGA-COAD'; 'STAD' = 'TCGA-STAD'; 'PRAD' = 'TCGA-PRAD'; 'THCA' = 'TCGA-THCA'; 'HNSC' = 'TCGA-HNSC' }
New-Item -ItemType Directory -Force -Path 'new\data\gdc_mirna' | Out-Null
foreach ($abbr in $cohorts.Keys) {
  $proj = $cohorts[$abbr]
  $filters = '{"op":"and","content":[' +
    '{"op":"in","content":{"field":"cases.project.project_id","value":["' + $proj + '"]}},' +
    '{"op":"in","content":{"field":"files.data_type","value":["miRNA Expression Quantification"]}},' +
    '{"op":"in","content":{"field":"files.experimental_strategy","value":["miRNA-Seq"]}},' +
    '{"op":"in","content":{"field":"cases.samples.sample_type","value":["Primary Tumor"]}}]}'
  $body = @{ filters = $filters; fields = 'file_id,file_name,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type'; size = '1000'; format = 'JSON' } | ConvertTo-Json
  $ok = $false
  for ($i = 0; $i -lt 4 -and -not $ok; $i++) {
    try {
      $r = Invoke-RestMethod -Uri 'https://api.gdc.cancer.gov/files' -Method Post -Body $body -ContentType 'application/json' -TimeoutSec 90
      $rows = foreach ($h in $r.data.hits) {
        $s = @($h.cases[0].samples) | Where-Object { $_.sample_type -eq 'Primary Tumor' } | Select-Object -First 1
        [pscustomobject]@{ file_id = $h.file_id; file_name = $h.file_name; case_id = $h.cases[0].submitter_id; sample_id = $s.submitter_id }
      }
      $rows | Export-Csv -Path "new\data\gdc_mirna\manifest_$abbr.csv" -NoTypeInformation -Encoding UTF8
      Write-Host "$abbr files=$($r.data.pagination.total) manifest_rows=$(@($rows).Count)"
      $ok = $true
    } catch { Write-Host "$abbr retry $i : $($_.Exception.Message)"; Start-Sleep -Seconds 4 }
  }
  if (-not $ok) { Write-Host "$abbr FAILED" }
}
