required_cran <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "dplyr",
  "readr",
  "tibble",
  "RColorBrewer"
)

required_bioc <- c(
  "DESeq2",
  "airway",
  "SummarizedExperiment",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "apeglm"
)

install_if_missing <- function() {
  missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_cran) > 0) {
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
}

install_if_missing()

suppressPackageStartupMessages({
  library(DESeq2)
  library(airway)
  library(SummarizedExperiment)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(readr)
  library(tibble)
  library(RColorBrewer)
})

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)

theme_publication <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey35"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

deg_alpha <- 0.05
deg_lfc <- 1

