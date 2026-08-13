vsd <- readRDS("data/processed/deseq2_vst.rds")
metadata <- as.data.frame(colData(vsd))
res_df <- read_csv("results/tables/deg_results.csv", show_col_types = FALSE)
metadata$sample_label <- paste0(ifelse(metadata$dex == "trt", "T", "U"), "-", metadata$cell)

pca_data <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_data$sample_label <- metadata[pca_data$name, "sample_label"]

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = dex, shape = cell, label = sample_label)) +
  geom_point(size = 3.2, alpha = 0.95) +
  geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c("untrt" = "#4C78A8", "trt" = "#F58518")) +
  labs(
    title = "PCA of airway RNA-seq samples",
    subtitle = "Variance-stabilized expression values",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Treatment",
    shape = "Cell line"
  ) +
  theme_publication()

ggsave("results/plots/pca.png", p_pca, width = 7, height = 5, dpi = 300)
ggsave("results/plots/pca.pdf", p_pca, width = 7, height = 5)

sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)
rownames(sample_dist_matrix) <- metadata$sample_label
colnames(sample_dist_matrix) <- rownames(sample_dist_matrix)

png("results/plots/sample_distance_heatmap.png", width = 1800, height = 1500, res = 250)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  main = "Sample distance heatmap"
)
dev.off()

pdf("results/plots/sample_distance_heatmap.pdf", width = 7, height = 6)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  main = "Sample distance heatmap"
)
dev.off()

top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 40) %>%
  pull(gene_id)

mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- metadata[colnames(mat), "sample_label"]
gene_labels <- res_df$gene_symbol[match(rownames(mat), res_df$gene_id)]
gene_labels[is.na(gene_labels) | gene_labels == ""] <- rownames(mat)[is.na(gene_labels) | gene_labels == ""]
rownames(mat) <- make.unique(gene_labels)

annotation_col <- metadata[, c("dex", "cell"), drop = FALSE]
rownames(annotation_col) <- metadata$sample_label

png("results/plots/top_deg_heatmap.png", width = 1800, height = 2100, res = 250)
pheatmap(
  mat,
  annotation_col = annotation_col,
  show_colnames = TRUE,
  fontsize_row = 7,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
  main = "Top 40 differentially expressed genes"
)
dev.off()

pdf("results/plots/top_deg_heatmap.pdf", width = 7, height = 8)
pheatmap(
  mat,
  annotation_col = annotation_col,
  show_colnames = TRUE,
  fontsize_row = 7,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
  main = "Top 40 differentially expressed genes"
)
dev.off()

volcano_df <- res_df %>%
  mutate(
    neg_log10_padj = -log10(padj),
    label = ifelse(significant & row_number() <= 12, coalesce(gene_symbol, gene_id), NA_character_)
  )

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, neg_log10_padj, color = direction)) +
  geom_point(alpha = 0.65, size = 1.35, na.rm = TRUE) +
  geom_vline(xintercept = c(-deg_lfc, deg_lfc), linetype = "dashed", color = "grey45") +
  geom_hline(yintercept = -log10(deg_alpha), linetype = "dashed", color = "grey45") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE, show.legend = FALSE) +
  scale_color_manual(values = c(
    "Up in treated" = "#D55E00",
    "Down in treated" = "#0072B2",
    "Not significant" = "grey70"
  )) +
  labs(
    title = "Volcano plot: dexamethasone treated vs untreated",
    subtitle = "Significance threshold: padj < 0.05 and |log2FC| >= 1",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value",
    color = "DEG status"
  ) +
  theme_publication()

ggsave("results/plots/volcano.png", p_volcano, width = 7, height = 5.5, dpi = 300)
ggsave("results/plots/volcano.pdf", p_volcano, width = 7, height = 5.5)

cat("Plots saved to results/plots\n")
