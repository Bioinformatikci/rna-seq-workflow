data("airway", package = "airway")

se <- airway
se$dex <- relevel(factor(se$dex), ref = "untrt")

counts <- assay(se)
metadata <- as.data.frame(colData(se)) %>%
  rownames_to_column("sample_id") %>%
  select(sample_id, cell, dex, albut, Run, avgLength, Experiment, Sample, BioSample)

write_csv(as.data.frame(counts) %>% rownames_to_column("gene_id"), "data/raw/airway_counts.csv")
write_csv(metadata, "data/raw/airway_metadata.csv")
saveRDS(se, "data/raw/airway_summarized_experiment.rds")

cat("Saved raw count matrix and sample metadata to data/raw\n")

