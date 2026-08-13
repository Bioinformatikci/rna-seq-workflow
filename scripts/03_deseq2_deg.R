qc <- readRDS("data/processed/qc_data.rds")

dds <- DESeqDataSetFromMatrix(
  countData = round(qc$counts),
  colData = qc$metadata,
  design = ~ cell + dex
)

dds <- DESeq(dds)

coef_name <- "dex_trt_vs_untrt"
res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    gene_symbol = mapIds(
      org.Hs.eg.db,
      keys = gene_id,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    ),
    significant = !is.na(padj) & padj < deg_alpha & abs(log2FoldChange) >= deg_lfc,
    direction = case_when(
      significant & log2FoldChange > 0 ~ "Up in treated",
      significant & log2FoldChange < 0 ~ "Down in treated",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(padj)

norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

vsd <- vst(dds, blind = FALSE)

write_csv(res_df, "results/tables/deg_results.csv")
write_csv(filter(res_df, significant), "results/tables/significant_deg_results.csv")
write_csv(norm_counts, "results/tables/normalized_counts.csv")
saveRDS(dds, "data/processed/deseq2_dds.rds")
saveRDS(vsd, "data/processed/deseq2_vst.rds")

capture.output(sessionInfo(), file = "results/session_info.txt")

cat("DESeq2 analysis complete. DEG tables saved to results/tables\n")

