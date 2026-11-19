#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(pheatmap)
    library(org.Hs.eg.db)
})

cat("=== SCRIPT 3: VISUALIZACIONES LOCALES Y REPARACIÓN DE ADN ===\n")

dir.create("plots", showWarnings = FALSE, recursive = TRUE)

dds <- readRDS("results/dds_fitted.rds")
res_df <- read.csv("results/deseq2_all_genes_results.csv")
degs <- read.csv("results/degs_significant_FDR0.05_FC1.5.csv")

norm_counts <- counts(dds, normalized = TRUE)
symbols <- mapIds(org.Hs.eg.db, keys = sub("\\..*", "", rownames(norm_counts)), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

valid_symbols <- ifelse(is.na(symbols), rownames(norm_counts), symbols)
rownames(norm_counts) <- make.unique(valid_symbols)

# 1. Heatmap Top DEGs con barra de anotación de muestras (Control vs Endometriosis)
top_genes <- degs %>% arrange(padj) %>% head(30) %>% pull(symbol)
top_genes <- intersect(top_genes, rownames(norm_counts))

if (length(top_genes) > 0) {
    mat_top <- norm_counts[top_genes, , drop = FALSE]
    mat_top_z <- t(scale(t(mat_top)))
    mat_top_z[is.na(mat_top_z)] <- 0

    # Crear data frame de anotación para las columnas
    df_annotation <- data.frame(
        Grupo = colData(dds)$condition,
        row.names = colnames(dds)
    )
    
    ann_colors <- list(
        Grupo = c("Control" = "#2b5c8f", "Endometriosis" = "#d95f02")
    )

    png("plots/05_heatmap_top_degs.png", width = 7.5, height = 8, units = "in", res = 300)
    pheatmap(mat_top_z, main = "Top 30 DEGs (Z-score)", cluster_cols = TRUE, show_colnames = TRUE, annotation_col = df_annotation, annotation_colors = ann_colors)
    invisible(dev.off())
}

# 2. Volcano Plot DDR con líneas de umbral y leyenda explícita
ddr_genes <- c("MSH6", "PARP2", "UNG", "TDG", "RAD51", "GEN1", "FAAP24", "PRKDC", "ATM", "TP53", "GADD45G", "DDIT3")
res_ddr <- res_df %>% filter(symbol %in% ddr_genes)

p_volcano_ddr <- ggplot() +
    geom_point(data = res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = "Fondo genómico"), alpha = 0.4, size = 1.2) +
    geom_point(data = res_ddr, aes(x = log2FoldChange, y = -log10(pvalue), color = "Genes DDR"), size = 3) +
    geom_text(data = res_ddr, aes(x = log2FoldChange, y = -log10(pvalue), label = symbol), vjust = -0.7, size = 3.5, fontface = "bold") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30") +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "gray30") +
    scale_color_manual(name = "Categoría", values = c("Fondo genómico" = "gray70", "Genes DDR" = "#d95f02")) +
    labs(title = "Volcano Plot: Genes de Reparación de ADN (DDR)", subtitle = "Líneas punteadas: FDR < 0.05 y |log2FC| ≥ 0.585", x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic"), legend.position = "right")

ggsave("plots/03b_volcano_ddr_genes.png", p_volcano_ddr, width = 7.5, height = 5, dpi = 300)

# 3. Barplots por Vías
pathways_genes <- list(
    "MMR"  = c("MLH1", "MSH3", "MSH6", "PCNA", "EXO1", "LIG1"),
    "BER"  = c("OGG1", "UNG", "LIG1", "FEN1", "POLB", "PNKP", "RFC1", "PCNA", "SMUG1", "MBD4", "TDG", "NEIL3", "PARP2", "APEX1"),
    "NER"  = c("DDB2", "RFC1", "PCNA", "GTF2H6", "LIG1", "FEN1", "CETN2", "RPA2", "RPA3", "GTF2H2", "MNAT1"),
    "HR"   = c("RAD51", "DMC1", "RAD54B", "RAD54L", "BLM", "BRCA2", "BARD1", "MRE11", "EXO1", "GEN1", "SLX1B", "RAD51B", "BRCA1", "EME2", "ATM"),
    "FA"   = c("FANCD2", "FANCI", "BRIP1", "FANCO", "FAAP20", "FAAP24", "BRCA2"),
    "NHEJ" = c("XRCC4", "DCLRE1C", "PNKP", "POLM", "POLL")
)

cond <- colData(dds)$condition
plot_data_list <- list()

for (path_name in names(pathways_genes)) {
    genes <- pathways_genes[[path_name]]
    present_genes <- intersect(genes, rownames(norm_counts))

    for (g in present_genes) {
        vals_ce <- norm_counts[g, cond == "Control"]
        vals_ee <- norm_counts[g, cond == "Endometriosis"]

        m_ce <- mean(vals_ce)
        s_ce <- if (length(vals_ce) > 1) sd(vals_ce) / sqrt(length(vals_ce)) else 0
        m_ee <- mean(vals_ee)
        s_ee <- if (length(vals_ee) > 1) sd(vals_ee) / sqrt(length(vals_ee)) else 0

        plot_data_list[[length(plot_data_list) + 1]] <- data.frame(Pathway = path_name, Gene = g, Condition = "CE", MeanExp = m_ce, SEM = s_ce)
        plot_data_list[[length(plot_data_list) + 1]] <- data.frame(Pathway = path_name, Gene = g, Condition = "EE", MeanExp = m_ee, SEM = s_ee)
    }
}

if (length(plot_data_list) > 0) {
    df_plot <- do.call(rbind, plot_data_list)
    df_plot$Condition <- factor(df_plot$Condition, levels = c("CE", "EE"))

    p_barplots <- ggplot(df_plot, aes(x = Gene, y = MeanExp, fill = Condition)) +
        geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
        geom_errorbar(aes(ymin = pmax(0, MeanExp - SEM), ymax = MeanExp + SEM), position = position_dodge(0.8), width = 0.25) +
        facet_wrap(~ Pathway, scales = "free", ncol = 3) +
        scale_fill_manual(values = c("CE" = "white", "EE" = "#d95f02"), labels = c("CE (Control)", "EE (Endometriosis)")) +
        labs(title = "Expresión de Genes de Reparación de ADN por Vías", x = "Gen", y = "Expresión Media Normalizada") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
              strip.background = element_rect(fill = "#2b5c8f"),
              strip.text = element_text(color = "white", face = "bold", size = 10),
              legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))

    ggsave("plots/04_fig4_dna_repair_barplots.png", p_barplots, width = 12, height = 8, dpi = 300)
}

# 4. Tabla Comparativa
key_genes_paper <- data.frame(
    symbol = c("MSH6", "PARP2", "UNG", "TDG", "GEN1", "RAD51", "FAAP24", "PRKDC"),
    via = c("MMR", "BER", "BER", "BER", "HR", "HR", "FA", "NHEJ/DDR"),
    p_value_paper = c("P = 0.047", "P < 0.05", "P = 0.03", "P = 0.04", "P = 0.04", "P = 0.04", "P = 0.02", "P < 0.05")
)

tabla_comp <- res_df %>%
    filter(symbol %in% key_genes_paper$symbol) %>%
    left_join(key_genes_paper, by = "symbol") %>%
    dplyr::select(symbol, via, log2FoldChange, pvalue, padj, p_value_paper) %>%
    arrange(pvalue)

write.csv(tabla_comp, "results/tabla_comparativa_genes_reparacion.csv", row.names = FALSE)

cat("Script 3 finalizado con éxito.\n")
