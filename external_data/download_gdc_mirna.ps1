# Batch-download GDC miRNA quantification files per cohort via POST /data {"ids":[...]}.
$ErrorActionPreference = 'Stop'
$cohorts = 'KIRC','LUAD','COAD','STAD','PRAD','THCA','HNSC'
$batchSize = 100
foreach ($abbr in $cohorts) {
  $manifest = Import-Csv "new\data\gdc_mirna\manifest_$abbr.csv"
  $rawDir = "new\data\gdc_mirna\raw_$abbr"
  New-Item -ItemType Directory -Force -Path $rawDir | Out-Null
  $have = @{}
  if (Test-Path "$rawDir\_downloaded.txt") { Get-Content "$rawDir\_downloaded.txt" | ForEach-Object { $have[$_] = $true } }
  $ids = @($manifest | Where-Object { -not $have.ContainsKey($_.file_id) } | ForEach-Object { $_.file_id })
  $total = $ids.Count
  Write-Host "$abbr to_download=$total"
  for ($start = 0; $start -lt $total; $start += $batchSize) {
    $chunk = @($ids[$start..([Math]::Min($start + $batchSize - 1, $total - 1))])
    $json = ConvertTo-Json -InputObject @{ ids = $chunk } -Compress
    $tmpTar = "$rawDir\_batch.tar"
    $done = $false
    for ($attempt = 1; $attempt -le 3 -and -not $done; $attempt++) {
      try {
        Invoke-WebRequest -Uri 'https://api.gdc.cancer.gov/data' -Method Post -Body $json -ContentType 'application/json' -OutFile $tmpTar -TimeoutSec 600
        tar -xf $tmpTar -C $rawDir 2>&1 | Out-Null
        $found = 0
        foreach ($id in $chunk) { if (Test-Path "$rawDir\$id") { $found++ } }
        if ($found -eq $chunk.Count) { $done = $true }
        else { Write-Host "$abbr batch@${start}: extracted $found/$($chunk.Count), attempt $attempt"; Start-Sleep -Seconds 5 }
      } catch { Write-Host "$abbr batch@${start} attempt ${attempt}: $($_.Exception.Message)"; Start-Sleep -Seconds 8 }
    }
    if (Test-Path $tmpTar) { Remove-Item $tmpTar -Force }
    if ($done) { Add-Content -Path "$rawDir\_downloaded.txt" -Value $chunk }
    else { Write-Host "$abbr batch@${start} FAILED - will retry at end"; Add-Content -Path "$rawDir\_failed_batches.txt" -Value $start }
  }
  $have2 = @(Get-ChildItem "$rawDir" -Directory -ErrorAction SilentlyContinue).Count
  Write-Host "$abbr extracted_dirs=$have2 expected=$($manifest.Count)"
}
Write-Host 'ALL COHORTS DONE'
