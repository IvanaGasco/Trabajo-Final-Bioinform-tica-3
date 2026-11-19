#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

cat("=== SCRIPT 1: ANÁLISIS EXPLORATORIO Y QC ===\n")

dir.create("plots", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

counts_file <- "../05_quantification/gene_counts_raw.txt"
if (!file.exists(counts_file)) {
  counts_file <- "../05_quantification/gene_counts.txt"
}
if (!file.exists(counts_file)) {
  stop("No se encontró la matriz de conteos en ../05_quantification/")
}

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
if (all(c("Chr", "Start", "End", "Strand", "Length") %in% colnames(counts))) {
  counts <- counts[, 6:ncol(counts)]
}

# Limpiar nombres de columnas
colnames(counts) <- gsub("\\.sorted\\.bam$", "", colnames(counts))
colnames(counts) <- gsub("\\.bam$", "", colnames(counts))
colnames(counts) <- gsub(".*/", "", colnames(counts))

samples <- colnames(counts)
condition <- ifelse(grepl("403|404|405|406", samples), "Control", "Endometriosis")

coldata <- data.frame(
  row.names = samples,
  condition = factor(condition, levels = c("Control", "Endometriosis"))
)
write.csv(coldata, "results/metadata.csv")

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

saveRDS(dds, "results/dds_fitted.rds")
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, "results/vsd_data.rds")

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.7, size = 3) +
  scale_color_manual(values = c("Control" = "#2b5c8f", "Endometriosis" = "#d95f02")) +
  labs(title = "PCA Samples (VST)", x = paste0("PC1: ", percentVar[1], "% var"), y = paste0("PC2: ", percentVar[2], "% var")) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("plots/01_pca_plot.png", p_pca, width = 7, height = 5, dpi = 300)

sampleDists <- dist(t(assay(vsd)))
png("plots/02_sample_distance_heatmap.png", width = 7, height = 6, units = "in", res = 300)
pheatmap(as.matrix(sampleDists),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")
dev.off()

cat("Script 1 finalizado con éxito.\n")
