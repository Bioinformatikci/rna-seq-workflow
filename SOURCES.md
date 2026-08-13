# Sources Checked

Internet sources checked on 2026-08-13:

- DESeq2 vignette: https://bioconductor.unipi.it/packages/3.19/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- Bioconductor RNA-seq gene-level workflow: https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
- Bioconductor airway package page: https://bioconductor.uib.no/packages/3.19/data/experiment/html/airway.html
- Airway workflow notes: https://bioconnector.github.io/workshops/r-rnaseq-airway.html

Key implementation choices:

- Use raw counts as DESeq2 input.
- Use `~ cell + dex` to estimate treatment effect while controlling for cell line.
- Use variance-stabilizing transformation for PCA and heatmaps.
- Report DEGs using `padj < 0.05` and `abs(log2FoldChange) >= 1`.

