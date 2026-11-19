#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(gprofiler2)
    library(stringr)
})

cat("=================================================================\n")
cat("=== PASO 8: ENRIQUECIMIENTO CON gprofiler2 (DISEÑO OPTIMIZADO) =\n")
cat("=================================================================\n\n")

dir.create("results/enrichment_lineaA", showWarnings = FALSE, recursive = TRUE)
dir.create("results/enrichment_lineaB", showWarnings = FALSE, recursive = TRUE)
dir.create("plots/enrichment_lineaA", showWarnings = FALSE, recursive = TRUE)
dir.create("plots/enrichment_lineaB", showWarnings = FALSE, recursive = TRUE)

file_a <- "../06_deseq2/results/deseq2_all_genes_results.csv"
file_b <- "../07_tuxedo_transcripts/results/t_test_fpkm_results.csv"

run_gprofiler <- function(genes, line_label, res_dir, plot_dir) {
    cat(sprintf("\n>>> [%s] Consultando g:Profiler API para %d genes...\n", line_label, length(genes)))
    
    gost_res <- tryCatch(
        gost(
            query = genes, 
            organism = "hsapiens", 
            ordered_query = FALSE, 
            multi_query = FALSE, 
            significant = TRUE, 
            exclude_iea = TRUE, 
            user_threshold = 0.05, 
            correction_method = "g_SCS", 
            sources = c("GO:BP", "KEGG", "REAC", "WP")
        ),
        error = function(e) {
            cat("Error de red / conexión con g:Profiler API:", e$message, "\n")
            return(NULL)
        }
    )
    
    if (!is.null(gost_res) && !is.null(gost_res$result) && nrow(gost_res$result) > 0) {
        df <- gost_res$result %>%
            select(any_of(c("source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "precision", "recall"))) %>%
            arrange(p_value)
            
        csv_out <- file.path(res_dir, "gprofiler_enrichment_results.csv")
        write.csv(df, csv_out, row.names = FALSE)
        cat(sprintf("Guardado CSV en: %s\n", csv_out))
        
        # Top 5 términos por fuente
        top_terms <- df %>%
            group_by(source) %>%
            slice_min(order_by = p_value, n = 5) %>%
            ungroup() %>%
            mutate(term_clean = stringr::str_wrap(term_name, width = 42))
            
        top_terms$term_clean <- factor(top_terms$term_clean, levels = rev(unique(top_terms$term_clean)))
        
        # Dotplot con Eje X ajustado y puntos aumentados
        p_dot <- ggplot(top_terms, aes(x = precision, y = term_clean, size = intersection_size, color = -log10(p_value))) +
            geom_point(alpha = 0.9) +
            facet_grid(source ~ ., scales = "free", space = "free_y") +
            scale_x_continuous(expand = expansion(mult = c(0.08, 0.15))) +
            scale_size_continuous(range = c(3.5, 8.5)) +
            scale_color_gradient(low = "#3182bd", high = "#de2d26") +
            labs(
                title = paste0("Vías Canónicas Enriquecidas - ", line_label),
                subtitle = "Enriquecimiento por g:Profiler (FDR < 0.05)",
                x = "Precisión (Genes coincidentes / Total en vía)",
                y = "Vía / Término Biológico",
                color = "-log10(p-value)",
                size = "N° Genes"
            ) +
            theme_bw(base_size = 11) +
            theme(
                plot.title = element_text(face = "bold", hjust = 0.5, size = 12.5),
                plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 10),
                axis.text.y = element_text(size = 8.5, face = "bold", color = "black"),
                axis.text.x = element_text(size = 9, color = "black"),
                strip.background = element_rect(fill = "#1f4e79"),
                strip.text = element_text(color = "white", face = "bold", size = 9.5),
                panel.spacing = unit(0.5, "lines")
            )
            
        plot_out <- file.path(plot_dir, "01_top_pathways_dotplot.png")
        ggsave(plot_out, p_dot, width = 10, height = 8.5, dpi = 300)
        cat(sprintf(" ✓ Gráfico optimizado guardado en: %s\n", plot_out))
        
        # Gráficos individuales por base de datos (con rango x ajustado)
        for (src in unique(top_terms$source)) {
            sub_df <- top_terms %>% filter(source == src)
            p_sub <- ggplot(sub_df, aes(x = precision, y = term_clean, size = intersection_size, color = -log10(p_value))) +
                geom_point(alpha = 0.9) +
                scale_x_continuous(expand = expansion(mult = c(0.08, 0.15))) +
                scale_size_continuous(range = c(4, 9)) +
                scale_color_gradient(low = "#3182bd", high = "#de2d26") +
                labs(
                    title = paste0("Vías Enriquecidas: ", src, " (", line_label, ")"),
                    x = "Precisión",
                    y = "",
                    color = "-log10(p-value)",
                    size = "N° Genes"
                ) +
                theme_bw(base_size = 11) +
                theme(
                    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
                    axis.text.y = element_text(size = 9, face = "bold", color = "black")
                )
                
            clean_src <- gsub(":", "_", src)
            ggsave(file.path(plot_dir, paste0("02_dotplot_", clean_src, ".png")), p_sub, width = 8.5, height = 4.8, dpi = 300)
        }
        cat(" ✓ Gráficos individuales guardados.\n")
    } else {
        cat(" - No se encontraron términos significativos.\n")
    }
}

# 1. LÍNEA A
if (file.exists(file_a)) {
    df_a <- read.csv(file_a)
    degs_a <- df_a %>%
        filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 0.585) %>%
        pull(symbol) %>%
        unique()
    run_gprofiler(degs_a, "Línea A (DESeq2 DEGs)", "results/enrichment_lineaA", "plots/enrichment_lineaA")
}

# 2. LÍNEA B
if (file.exists(file_b)) {
    df_b <- read.csv(file_b)
    lfc_col <- if ("log2FoldChange" %in% colnames(df_b)) "log2FoldChange" else "log2FC"
    degs_b <- df_b %>%
        filter((!is.na(p_value) & p_value < 0.05) | abs(.data[[lfc_col]]) >= 0.585) %>%
        pull(symbol) %>%
        unique()
    run_gprofiler(degs_b, "Línea B (Tuxedo FPKM)", "results/enrichment_lineaB", "plots/enrichment_lineaB")
}

cat("\n=== PROCESAMIENTO FINALIZADO CON EXITO ===\n")
