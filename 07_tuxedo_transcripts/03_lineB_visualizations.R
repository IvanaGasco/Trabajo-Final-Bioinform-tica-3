#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    library(scales)
})

cat("=== LÍNEA B - SCRIPT 3: VISUALIZACIONES COMPLETAS (6 VÍAS + GENES CLAVE + TOP 10) ===\n")

dir.create("plots", showWarnings = FALSE, recursive = TRUE)

fpkm_file <- "results/isoforms_fpkm_repair_genes.csv"
ttest_file <- "results/t_test_fpkm_results.csv"

if (!file.exists(fpkm_file)) {
    stop("ERROR: No se encontró 'results/isoforms_fpkm_repair_genes.csv'. Ejecutá los pasos previos.")
}

fpkm_df <- read.csv(fpkm_file)

# Vías de reparación de ADN y sus genes
pathways_genes <- list(
    "MMR"  = c("MLH1", "MSH3", "MSH6", "PCNA", "EXO1", "LIG1"),
    "BER"  = c("OGG1", "UNG", "LIG1", "FEN1", "POLB", "PNKP", "RFC1", "PCNA", "SMUG1", "MBD4", "TDG", "NEIL3", "PARP2", "APEX1"),
    "NER"  = c("DDB2", "RFC1", "PCNA", "GTF2H6", "LIG1", "FEN1", "CETN2", "RPA2", "RPA3", "GTF2H2", "MNAT1"),
    "HR"   = c("RAD51", "DMC1", "RAD54B", "RAD54L", "BLM", "BRCA2", "BARD1", "MRE11", "EXO1", "GEN1", "SLX1B", "RAD51B", "BRCA1", "EME2", "ATM"),
    "FA"   = c("FANCD2", "FANCI", "BRIP1", "FANCO", "FAAP20", "FAAP24", "BRCA2"),
    "NHEJ" = c("XRCC4", "DCLRE1C", "PNKP", "POLM", "POLL")
)

pathway_df <- do.call(rbind, lapply(names(pathways_genes), function(p) {
    data.frame(gene_symbol = pathways_genes[[p]], Pathway = p, stringsAsFactors = FALSE)
}))

fpkm_pathways <- fpkm_df %>% inner_join(pathway_df, by = "gene_symbol")

summary_all <- fpkm_pathways %>%
    group_by(Pathway, gene_symbol, condition) %>%
    summarise(
        mean_FPKM = mean(FPKM, na.rm = TRUE),
        sem_FPKM = ifelse(n() > 1, sd(FPKM, na.rm = TRUE) / sqrt(n()), 0),
        .groups = "drop"
    )

summary_all$condition <- factor(summary_all$condition, levels = c("Control", "Endometriosis"))

# ------------------------------------------------------------------------------
# Transformación Asimétrica por Tramos (Piecewise Transformation)
# Ocupa del 0 al 50% de la altura vertical para el rango FPKM 0 a 2.
# ------------------------------------------------------------------------------
piecewise_0_2 <- trans_new(
    name = "piecewise_0_2",
    transform = function(x) {
        x_clean <- pmax(0, x)
        ifelse(x_clean <= 2, x_clean / 4, 0.5 + log2(pmax(1, x_clean / 2)) * 0.15)
    },
    inverse = function(y) {
        y_clean <- pmax(0, y)
        ifelse(y_clean <= 0.5, y_clean * 4, 2 * (2^((y_clean - 0.5) / 0.15)))
    }
)

# ------------------------------------------------------------------------------
# 1. Gráfico Completo de las 6 Vías (Eje Y Asimétrico + Etiquetas en todos los paneles)
# ------------------------------------------------------------------------------
p_all_asym <- ggplot(summary_all, aes(x = gene_symbol, y = mean_FPKM, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
    geom_errorbar(aes(ymin = pmax(0, mean_FPKM - sem_FPKM), ymax = mean_FPKM + sem_FPKM), position = position_dodge(0.8), width = 0.25) +
    facet_wrap(~ Pathway, scales = "free", ncol = 3) +
    scale_y_continuous(
        trans = piecewise_0_2,
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.1)),
        breaks = c(0, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50, 100)
    ) +
    scale_fill_manual(values = c("Control" = "white", "Endometriosis" = "#d95f02"), labels = c("Control (CE)", "Endometriosis (EE)")) +
    labs(
        title = "Línea B: Expresión de Vías de Reparación de ADN (FPKM)",
        subtitle = "Eje Y asimétrico: Rango 0-2 FPKM expandido al 50% de altura para visibilizar valores bajos",
        x = "Gen",
        y = "FPKM Medio (Escala Asimétrica 0-2)"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic"),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5, face = "bold"),
        strip.background = element_rect(fill = "#2b5c8f"),
        strip.text = element_text(color = "white", face = "bold", size = 10)
    )

ggsave("plots/02_lineB_fpkm_all_pathways.png", p_all_asym, width = 13, height = 8.5, dpi = 300)
cat("✓ 1. Gráfico completo de 6 vías con eje Y asimétrico guardado en 'plots/02_lineB_fpkm_all_pathways.png'\n")

# ------------------------------------------------------------------------------
# 2. Gráfico Acotado: Genes Clave Destacados del Paper
# ------------------------------------------------------------------------------
key_genes_paper <- c("MSH6", "PARP2", "UNG", "TDG", "RAD51", "GEN1", "FAAP24", "PRKDC")
summary_key <- summary_all %>% filter(gene_symbol %in% key_genes_paper)

p_key <- ggplot(summary_key, aes(x = gene_symbol, y = mean_FPKM, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
    geom_errorbar(aes(ymin = pmax(0, mean_FPKM - sem_FPKM), ymax = mean_FPKM + sem_FPKM), position = position_dodge(0.8), width = 0.25) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    scale_fill_manual(values = c("Control" = "white", "Endometriosis" = "#d95f02"), labels = c("Control (CE)", "Endometriosis (EE)")) +
    labs(
        title = "Línea B: Expresión de Genes Clave de Reparación (FPKM)",
        subtitle = "Medias de FPKM (± SEM) para genes clave representativos",
        x = "Gen",
        y = "FPKM Medio"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
    )

ggsave("plots/01_lineB_fpkm_barplots.png", p_key, width = 8, height = 5, dpi = 300)
cat("✓ 2. Gráfico de genes clave guardado en 'plots/01_lineB_fpkm_barplots.png'\n")

# ------------------------------------------------------------------------------
# 3. Gráfico Top 10 Genes Más Significativos de la Prueba t
# ------------------------------------------------------------------------------
if (file.exists(ttest_file)) {
    ttest_res <- read.csv(ttest_file)
    top10_genes <- ttest_res %>% filter(!is.na(p_value)) %>% arrange(p_value) %>% head(10) %>% pull(symbol)
    
    summary_top10 <- summary_all %>% filter(gene_symbol %in% top10_genes)
    
    p_top10 <- ggplot(summary_top10, aes(x = reorder(gene_symbol, mean_FPKM, max), y = mean_FPKM, fill = condition)) +
        geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
        geom_errorbar(aes(ymin = pmax(0, mean_FPKM - sem_FPKM), ymax = mean_FPKM + sem_FPKM), position = position_dodge(0.8), width = 0.25) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
        scale_fill_manual(values = c("Control" = "white", "Endometriosis" = "#d95f02"), labels = c("Control (CE)", "Endometriosis (EE)")) +
        labs(
            title = "Línea B: Top 10 Genes con Menor p-valor (FPKM)",
            subtitle = "Medias de FPKM (± SEM) ordenados por menor p-valor de la prueba t",
            x = "Gen",
            y = "FPKM Medio"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "top",
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
        )
        
    ggsave("plots/03_lineB_top10_pvalue_barplots.png", p_top10, width = 8, height = 5, dpi = 300)
    cat("✓ 3. Gráfico Top 10 por p-valor guardado en 'plots/03_lineB_top10_pvalue_barplots.png'\n")
}

cat("=== TODAS LAS VISUALIZACIONES DEL SCRIPT 3 COMPLETADAS EXITOSAMENTE ===\n")
