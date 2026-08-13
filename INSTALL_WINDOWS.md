# Windows Install Notes

This project needs R because the core differential expression method is DESeq2 from Bioconductor.

## Check Environment

From the project root:

```powershell
.\scripts\check_environment.ps1
```

## Install R

If `Rscript` is not found, install R:

```powershell
winget install RProject.R
```

Then close and reopen PowerShell or Codex so the PATH refreshes.

## Run Analysis

```powershell
Rscript run_analysis.R
```

The first run may take several minutes because R packages are installed automatically.

