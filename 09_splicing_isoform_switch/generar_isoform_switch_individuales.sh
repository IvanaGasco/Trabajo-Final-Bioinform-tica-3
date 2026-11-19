#!/bin/bash
# ==============================================================================
# Script: generar_isoform_switch_individuales.sh (Paso 9)
# Descripción: Calcula la Fracción de Isoforma (IF %) y genera Stacked Barplots
# individuales para PAQR5, MT1F, MT1X, S100P, GADD45G, DDIT3, etc.
# ==============================================================================

mkdir -p plots
mkdir -p results

cat << 'R_EOF' | Rscript -
suppressPackageStartupMessages({
    if (!require("dplyr", quietly = TRUE)) install.packages("dplyr", repos="https://cloud.r-project.org")
    if (!require("tidyr", quietly = TRUE)) install.packages("tidyr", repos="https://cloud.r-project.org")
    if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos="https://cloud.r-project.org")
    if (!require("scales", quietly = TRUE)) install.packages("scales", repos="https://cloud.r-project.org")
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(scales)
})

cat("=================================================================\n")
cat("=== GENERACIÓN DE ISOFORM SWITCH INDIVIDUAL (PANEL AMPLIADO) ===\n")
cat("=================================================================\n\n")

# 1. Rutas posibles de los datos de FPKM
file_ampliado <- "../07_tuxedo_transcripts/results/isoforms_fpkm_all_genes.csv"
file_repair <- "../07_tuxedo_transcripts/results/isoforms_fpkm_repair_genes.csv"
file_local <- "results/isoform_fraction_all_transcripts.csv"

if (file.exists(file_ampliado)) {
    cat("➜ Cargando FPKMs ampliados desde:", file_ampliado, "\n")
    df_raw <- read.csv(file_ampliado)
} else if (file.exists(file_repair)) {
    cat("➜ Cargando FPKMs desde:", file_repair, "\n")
    df_raw <- read.csv(file_repair)
} else {
    stop("Error: No se encontró ningún archivo de FPKM en el Paso 7.")
}

# 2. Promediar FPKM por transcrito y condición
df_summary <- df_raw %>%
    group_by(gene_symbol, transcript_id, condition) %>%
    summarise(mean_fpkm = mean(as.numeric(FPKM), na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = mean_fpkm, values_fill = 0)

cond_cols <- setdiff(colnames(df_summary), c("gene_symbol", "transcript_id"))
col_ce <- cond_cols[grep("Control|CE", cond_cols, ignore.case = TRUE)]
col_ee <- cond_cols[grep("Endometriosis|EE|Case", cond_cols, ignore.case = TRUE)]

df_wide <- df_summary %>%
    rename(
        symbol = gene_symbol,
        mean_CE = all_of(col_ce),
        mean_EE = all_of(col_ee)
    )

# 3. Expresión Total del Gen y Fracción de Isoforma (IF)
df_totals <- df_wide %>%
    group_by(symbol) %>%
    summarise(
        total_fpkm_CE = sum(mean_CE, na.rm = TRUE),
        total_fpkm_EE = sum(mean_EE, na.rm = TRUE),
        n_isoforms = n(),
        .groups = "drop"
    )

df_switch <- df_wide %>%
    inner_join(df_totals, by = "symbol") %>%
    filter(total_fpkm_CE > 0.0001 | total_fpkm_EE > 0.0001) %>%
    mutate(
        IF_CE = ifelse(total_fpkm_CE > 0, mean_CE / total_fpkm_CE, 0),
        IF_EE = ifelse(total_fpkm_EE > 0, mean_EE / total_fpkm_EE, 0),
        dIF = IF_EE - IF_CE
    )

# Guardar la tabla recalculada
write.csv(df_switch, "results/isoform_fraction_all_transcripts.csv", row.names = FALSE)

# 4. Lista de genes objetivo a graficar
genes_objetivo <- c("PAQR5", "MT1F", "MT1X", "S100P", "GADD45G", "DDIT3", "BRCA1", "BLM", "RAD51", "RAD54L", "MSH3")
genes_presentes <- intersect(genes_objetivo, unique(df_switch$symbol))

cat(paste0("➜ Genes encontrados y listos para graficar (", length(genes_presentes), "): ", paste(genes_presentes, collapse = ", "), "\n\n"))

# 5. Generación de gráficos individuales
for (gene in genes_presentes) {
    df_gene <- df_switch %>% filter(symbol == gene)
    
    if (nrow(df_gene) <= 1) {
        cat(paste0(gene, " posee solo 1 isoforma expresada. Se omite el gráfico de proporciones.\n"))
        next
    }
    
    # Agrupar isoformas con contribución < 5% como "Otras (<5%)"
    df_gene_prep <- df_gene %>%
        mutate(
            max_IF = pmax(IF_CE, IF_EE),
            tx_label = ifelse(max_IF < 0.05, "Otras isoformas (<5%)", transcript_id)
        ) %>%
        group_by(symbol, tx_label) %>%
        summarise(
            IF_CE = sum(IF_CE),
            IF_EE = sum(IF_EE),
            .groups = "drop"
        ) %>%
        pivot_longer(cols = c(IF_CE, IF_EE), names_to = "Condition", values_to = "IF") %>%
        mutate(Condition = factor(ifelse(Condition == "IF_CE", "Control", "Endometriosis"), levels = c("Control", "Endometriosis")))
        
    p_gene <- ggplot(df_gene_prep, aes(x = Condition, y = IF, fill = tx_label)) +
        geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.4, width = 0.55) +
        scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
        scale_fill_brewer(palette = "Set3") +
        labs(
            title = paste0("Redistribución de Isoformas: ", gene),
            subtitle = paste0("Fracción de Isoforma (IF %) | Control (CE-MS) vs. Endometriosis (EE-MS)"),
            x = "Grupo Biológico",
            y = "Proporción de Expresión (%)",
            fill = "Transcrito Ensembl"
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 10, color = "gray30"),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(face = "bold", color = "black", size = 11),
            legend.position = "right",
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 9),
            panel.grid.major.x = element_blank()
        )
        
    outfile <- paste0("plots/isoform_switch_", gene, ".png")
    ggsave(outfile, p_gene, width = 6.5, height = 5, dpi = 300)
    cat(paste0(" ✓ Guardado exitosamente: ", outfile, "\n"))
}

cat("\n=================================================================\n")
cat("=== PROCESO COMPLETADO EXITOSAMENTE ===\n")
cat("=================================================================\n")
R_EOF
