#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(org.Hs.eg.db)
})

cat("=== SCRIPT 2: EXPRESIÓN DIFERENCIAL (DESeq2) ===\n")

dds <- readRDS("results/dds_fitted.rds")
res <- results(dds, contrast = c("condition", "Endometriosis", "Control"))

res_shrunk <- tryCatch({
    lfcShrink(dds, coef = "condition_Endometriosis_vs_Control", type = "normal")
}, error = function(e) {
    res
})

res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)

symbols <- mapIds(org.Hs.eg.db, keys = sub("\\..*", "", res_df$gene_id), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res_df$symbol <- ifelse(is.na(symbols), res_df$gene_id, symbols)

write.csv(res_df, "results/deseq2_all_genes_results.csv", row.names = FALSE)

degs <- res_df %>% filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 0.585)
write.csv(degs, "results/degs_significant_FDR0.05_FC1.5.csv", row.names = FALSE)

res_df <- res_df %>% mutate(diff_expressed = case_when(
    padj < 0.05 & log2FoldChange >= 0.585 ~ "UP",
    padj < 0.05 & log2FoldChange <= -0.585 ~ "DOWN",
    TRUE ~ "NO"
))

# Volcano Plot General
p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = diff_expressed)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30") +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "gray30") +
    scale_color_manual(values = c("DOWN" = "#3182bd", "NO" = "gray70", "UP" = "#de2d26")) +
    labs(title = "Volcano Plot: Endometriosis vs Control", x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("plots/03_volcano_plot.png", p_volcano, width = 7, height = 5, dpi = 300)

# MA Plot acotado de -2.5 a 2.5
png("plots/04_ma_plot.png", width = 7, height = 5, units = "in", res = 300)
plotMA(res_shrunk, main = "MA Plot: DESeq2 (LFC Shrinkage)", ylim = c(-2.5, 2.5))
dev.off()

cat("Script 2 finalizado con éxito.\n")
