# RNA-seq Differential Expression Summary

## Dataset

This analysis uses the Bioconductor airway dataset from Himes et al. 2014. The experiment profiles human airway smooth muscle cells treated with dexamethasone and matched untreated controls.

## Analysis Design

- Method: DESeq2 negative binomial model
- Design formula: `~ cell + dex`
- Contrast: dexamethasone treated vs untreated
- DEG threshold: `padj < 0.05` and `abs(log2FoldChange) >= 1`
- Transformation for PCA and heatmaps: variance-stabilizing transformation

## QC Summary

- Samples analyzed: 8
- Genes before filtering: 63677
- Genes after filtering: 16139
- Filtering rule: keep genes with count >= 10 in at least 4 samples

## DEG Summary

- Significant DEGs: 749
- Up in treated: 411
- Down in treated: 338

## Top Significant Genes

| gene_id | gene_symbol | log2FoldChange | padj | direction |
| --- | --- | ---: | ---: | --- |
| ENSG00000152583 | SPARCL1 | 4.555676 | 5.530127e-132 | Up in treated |
| ENSG00000165995 | CACNB2 | 3.276456 | 5.530127e-132 | Up in treated |
| ENSG00000120129 | DUSP1 | 2.935137 | 1.248569e-124 | Up in treated |
| ENSG00000101347 | SAMHD1 | 3.754493 | 2.909302e-124 | Up in treated |
| ENSG00000189221 | MAOA | 3.337353 | 3.538038e-119 | Up in treated |
| ENSG00000211445 | GPX3 | 3.711558 | 1.344901e-106 | Up in treated |
| ENSG00000157214 | STEAP2 | 1.963647 | 3.955211e-101 | Up in treated |
| ENSG00000162614 | NEXN | 2.022927 |  1.212676e-96 | Up in treated |
| ENSG00000125148 | MT2A | 2.196818 |  1.961381e-91 | Up in treated |
| ENSG00000154734 | ADAMTS1 | 2.330306 |  4.017702e-85 | Up in treated |

## Deliverables

- `results/tables/deg_results.csv`: full differential expression table
- `results/tables/significant_deg_results.csv`: filtered significant DEG table
- `results/tables/normalized_counts.csv`: DESeq2-normalized counts
- `results/plots/pca.png`: PCA plot for sample-level QC
- `results/plots/sample_distance_heatmap.png`: sample distance heatmap
- `results/plots/top_deg_heatmap.png`: heatmap of top differentially expressed genes
- `results/plots/volcano.png`: volcano plot for DEG overview

## Interpretation

The workflow checks sample-level structure with PCA and distance-based clustering before interpreting differential expression. The model adjusts for donor/cell-line effects while estimating the treatment effect of dexamethasone. Genes passing the stated FDR and effect-size thresholds are reported as candidate treatment-responsive genes. These results are suitable for follow-up pathway analysis, target prioritization, or validation planning, but should not be over-interpreted without experimental context.
