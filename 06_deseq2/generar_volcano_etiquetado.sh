#!/bin/bash
# ==============================================================================
# Script: generar_volcano_etiquetado.sh
# Descripción: Genera el Volcano Plot global (DESeq2 - Línea A) con etiquetas
# automáticas (ggrepel) para los DEGs significativos y genes DDR clave.
# ==============================================================================

# Crear carpeta de salidas si no existe
mkdir -p plots

cat << 'R_EOF' | Rscript -
suppressPackageStartupMessages({
    if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos="https://cloud.r-project.org")
    if (!require("dplyr", quietly = TRUE)) install.packages("dplyr", repos="https://cloud.r-project.org")
    if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel", repos="https://cloud.r-project.org")
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
})

cat("=== Cargando datos y generando Volcano Plot Etiquetado ===\n")

# 1. Leer resultados completos de DESeq2
if (!file.exists("results/deseq2_all_genes_results.csv")) {
    stop("Error: No se encontró el archivo results/deseq2_all_genes_results.csv")
}
res_df <- read.csv("results/deseq2_all_genes_results.csv")

# 2. Definir estado de expresión (FDR < 0.05 y |log2FC| >= 0.585)
res_df <- res_df %>%
    mutate(diff_expressed = case_when(
        padj < 0.05 & log2FoldChange >= 0.585 ~ "UP",
        padj < 0.05 & log2FoldChange <= -0.585 ~ "DOWN",
        TRUE ~ "NO"
    ))

# 3. Seleccionar genes para etiquetar:
# - Todos los DEGs significativos
# - Genes DDR / respuesta a estrés clave (GADD45G, DDIT3, MT1X)
genes_a_etiquetar <- res_df %>%
    filter(diff_expressed %in% c("UP", "DOWN") | symbol %in% c("GADD45G", "DDIT3", "MT1X"))

# 4. Construir el gráfico
p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = diff_expressed)) +
    # Puntos de fondo
    geom_point(alpha = 0.55, size = 1.5) +
    # Umbrales estadísticos (p = 0.05 y log2FC = 0.585)
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "gray40", linewidth = 0.6) +
    # Esquema de colores
    scale_color_manual(values = c("DOWN" = "#3182bd", "NO" = "gray75", "UP" = "#de2d26")) +
    # Repel de texto sin solapamiento
    geom_text_repel(
        data = genes_a_etiquetar,
        aes(label = symbol),
        size = 3.6,
        fontface = "bold",
        box.padding = 0.4,
        point.padding = 0.3,
        segment.color = "gray30",
        segment.size = 0.4,
        max.overlaps = Inf,
        show.legend = FALSE
    ) +
    # Títulos y etiquetas
    labs(
        title = "Volcano Plot: Endometriosis vs Control (Línea A)",
        x = expression(log[2] ~ "Fold Change"),
        y = expression(-log[10] ~ "(" * p * "-value)"),
        color = "Expresión"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title = element_text(face = "bold"),
        legend.position = "right"
    )

# 5. Guardar imagen final en plots/
output_path <- "plots/03_volcano_plot_etiquetado.png"
ggsave(output_path, p_volcano, width = 8, height = 6, dpi = 300)
cat(paste0("¡Proceso finalizado exitosamente! Imagen guardada en: ", output_path, "\n"))
R_EOF
