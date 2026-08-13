$ErrorActionPreference = "Stop"

Write-Host "Checking RNA-seq workflow environment..."

$rscript = Get-Command Rscript -ErrorAction SilentlyContinue

if (-not $rscript) {
  Write-Host "Rscript was not found in PATH."
  Write-Host ""
  Write-Host "Install R, then reopen PowerShell/Codex and run:"
  Write-Host "  Rscript run_analysis.R"
  Write-Host ""
  Write-Host "Windows option:"
  Write-Host "  winget install RProject.R"
  exit 1
}

Write-Host "Found Rscript: $($rscript.Source)"
& Rscript --version

Write-Host ""
Write-Host "Environment check passed."

