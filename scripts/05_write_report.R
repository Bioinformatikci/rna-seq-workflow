qc <- read_csv("results/tables/qc_summary.csv", show_col_types = FALSE)
res_df <- read_csv("results/tables/deg_results.csv", show_col_types = FALSE)

samples <- qc$value[qc$metric == "samples"]
genes_before <- qc$value[qc$metric == "genes_before_filtering"]
genes_after <- qc$value[qc$metric == "genes_after_filtering"]
sig_total <- sum(res_df$significant, na.rm = TRUE)
sig_up <- sum(res_df$direction == "Up in treated", na.rm = TRUE)
sig_down <- sum(res_df$direction == "Down in treated", na.rm = TRUE)

top_table <- res_df %>%
  filter(significant) %>%
  select(gene_id, gene_symbol, log2FoldChange, padj, direction) %>%
  slice_head(n = 10)

top_lines <- apply(top_table, 1, function(x) {
  paste0("| ", paste(x, collapse = " | "), " |")
})

report <- c(
  "# RNA-seq Differential Expression Summary",
  "",
  "## Dataset",
  "",
  "This analysis uses the Bioconductor airway dataset from Himes et al. 2014. The experiment profiles human airway smooth muscle cells treated with dexamethasone and matched untreated controls.",
  "",
  "## Analysis Design",
  "",
  "- Method: DESeq2 negative binomial model",
  "- Design formula: `~ cell + dex`",
  "- Contrast: dexamethasone treated vs untreated",
  "- DEG threshold: `padj < 0.05` and `abs(log2FoldChange) >= 1`",
  "- Transformation for PCA and heatmaps: variance-stabilizing transformation",
  "",
  "## QC Summary",
  "",
  paste0("- Samples analyzed: ", samples),
  paste0("- Genes before filtering: ", genes_before),
  paste0("- Genes after filtering: ", genes_after),
  "- Filtering rule: keep genes with count >= 10 in at least 4 samples",
  "",
  "## DEG Summary",
  "",
  paste0("- Significant DEGs: ", sig_total),
  paste0("- Up in treated: ", sig_up),
  paste0("- Down in treated: ", sig_down),
  "",
  "## Top Significant Genes",
  "",
  "| gene_id | gene_symbol | log2FoldChange | padj | direction |",
  "| --- | --- | ---: | ---: | --- |",
  top_lines,
  "",
  "## Deliverables",
  "",
  "- `results/tables/deg_results.csv`: full differential expression table",
  "- `results/tables/significant_deg_results.csv`: filtered significant DEG table",
  "- `results/tables/normalized_counts.csv`: DESeq2-normalized counts",
  "- `results/plots/pca.png`: PCA plot for sample-level QC",
  "- `results/plots/sample_distance_heatmap.png`: sample distance heatmap",
  "- `results/plots/top_deg_heatmap.png`: heatmap of top differentially expressed genes",
  "- `results/plots/volcano.png`: volcano plot for DEG overview",
  "",
  "## Interpretation",
  "",
  "The workflow checks sample-level structure with PCA and distance-based clustering before interpreting differential expression. The model adjusts for donor/cell-line effects while estimating the treatment effect of dexamethasone. Genes passing the stated FDR and effect-size thresholds are reported as candidate treatment-responsive genes. These results are suitable for follow-up pathway analysis, target prioritization, or validation planning, but should not be over-interpreted without experimental context."
)

writeLines(report, "reports/summary_report.md")

cat("Summary report saved to reports/summary_report.md\n")

