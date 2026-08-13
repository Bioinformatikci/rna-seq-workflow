counts <- read_csv("data/raw/airway_counts.csv", show_col_types = FALSE) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

metadata <- read_csv("data/raw/airway_metadata.csv", show_col_types = FALSE) %>%
  as.data.frame()
rownames(metadata) <- metadata$sample_id

counts <- counts[, rownames(metadata)]
metadata$dex <- relevel(factor(metadata$dex), ref = "untrt")
metadata$cell <- factor(metadata$cell)

keep <- rowSums(counts >= 10) >= 4
filtered_counts <- counts[keep, ]

qc_summary <- tibble(
  metric = c(
    "samples",
    "genes_before_filtering",
    "genes_after_filtering",
    "filter_rule",
    "contrast"
  ),
  value = c(
    ncol(counts),
    nrow(counts),
    nrow(filtered_counts),
    "keep genes with count >= 10 in at least 4 samples",
    "treated vs untreated dexamethasone"
  )
)

write_csv(qc_summary, "results/tables/qc_summary.csv")
write_csv(as.data.frame(filtered_counts) %>% rownames_to_column("gene_id"), "data/processed/airway_counts_filtered.csv")
saveRDS(list(counts = filtered_counts, metadata = metadata), "data/processed/qc_data.rds")

cat("QC complete. Filtered count matrix saved to data/processed\n")

