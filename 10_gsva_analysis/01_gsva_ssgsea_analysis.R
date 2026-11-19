#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(GSVA)
    library(pheatmap)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
})

cat("=================================================================\n")
cat("=== ANALISIS ADICIONAL 2: GSVA & ssGSEA PATHWAY SCORES ==========\n")
cat("=================================================================\n\n")

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

# 1. Cargamos la matriz VST de DESeq2 (Paso 6)
vsd_file <- "../06_deseq2/results/vsd_data.rds"

if (!file.exists(vsd_file)) {
    stop("No se encontró la matriz VST en: ", vsd_file)
}

cat("✓ Cargando matriz VST desde:", vsd_file, "\n")
vsd <- readRDS(vsd_file)
expr_matrix <- assay(vsd)
metadata <- as.data.frame(colData(vsd))

cat(sprintf("✓ Matriz cargada: %d genes x %d muestras\n", nrow(expr_matrix), ncol(expr_matrix)))

# 2. Mapeo local de IDs Ensembl a Gene Symbols usando org.Hs.eg.db
cat("ℹ Mapeando IDs Ensembl a Gene Symbols mediante org.Hs.eg.db...\n")
clean_ensg <- sub("\\..*", "", rownames(expr_matrix))

map_df <- tryCatch({
    AnnotationDbi::select(org.Hs.eg.db, keys = clean_ensg, keytype = "ENSEMBL", columns = c("SYMBOL"))
}, error = function(e) {
    NULL
})

if (!is.null(map_df) && nrow(map_df) > 0) {
    map_df <- map_df[!is.na(map_df$SYMBOL) & map_df$SYMBOL != "", ]
    map_df <- map_df[!duplicated(map_df$ENSEMBL), ]
    
    symbol_map <- setNames(map_df$SYMBOL, map_df$ENSEMBL)
    mapped_symbols <- symbol_map[clean_ensg]
    mapped_symbols <- ifelse(is.na(mapped_symbols) | mapped_symbols == "", clean_ensg, mapped_symbols)
} else {
    mapped_symbols <- clean_ensg
}

# Consolidación de duplicados seleccionando el de mayor expresión promedio
df_temp <- as.data.frame(expr_matrix)
df_temp$Symbol <- mapped_symbols
df_temp$MeanExpr <- rowMeans(expr_matrix)

df_clean <- df_temp %>%
    group_by(Symbol) %>%
    arrange(desc(MeanExpr)) %>%
    slice(1) %>%
    ungroup()

# Seleccionamos únicamente las columnas de muestras de forma explícita sin conflicto de librerías
sample_cols <- setdiff(colnames(df_clean), c("Symbol", "MeanExpr"))
expr_matrix <- as.matrix(df_clean[, sample_cols])
rownames(expr_matrix) <- df_clean$Symbol

cat(sprintf("Matriz convertida a Símbolos: %d genes únicos x %d muestras\n\n", nrow(expr_matrix), ncol(expr_matrix)))

# 3. Definición de firmas DDR
ddr_gene_sets <- list(
    Homologous_Recombination = c("RAD51", "BRCA1", "BRCA2", "RAD54L", "EME2", "DMC1", "ATM", "BARD1", "MRE11", "RAD50", "NBN"),
    Base_Excision_Repair = c("APEX1", "PARP2", "TDG", "POLB", "OGG1", "LIG1", "FEN1", "MBD4", "PNKP"),
    Nucleotide_Excision_Repair = c("CETN2", "DDB2", "XPC", "ERCC2", "PCNA", "RFC1", "RPA1", "RPA2"),
    Mismatch_Repair = c("MSH6", "MSH3", "MLH1", "PMS2", "EXO1", "PCNA"),
    Non_Homologous_End_Joining = c("PRKDC", "XRCC5", "XRCC6", "LIG4", "DCLRE1C", "POLL", "POLM"),
    Fanconi_Anemia = c("FANCD2", "FANCI", "FAAP20", "FAAP24", "BRIP1", "FANCA", "FANCC"),
    p53_GADD45_Signaling = c("TP53", "GADD45G", "CDKN1A", "DDIT3", "MDM2", "CHEK2")
)

# 4. Cálculo de GSVA Scores
cat("Calculando GSVA Scores por muestra individual...\n")
param <- gsvaParam(exprData = expr_matrix, geneSets = ddr_gene_sets)
res_obj <- gsva(param)

if (inherits(res_obj, "SummarizedExperiment")) {
    gsva_result <- assay(res_obj)
} else {
    gsva_result <- as.matrix(res_obj)
}

write.csv(gsva_result, "results/gsva_pathway_scores_by_sample.csv")
cat("Matriz de GSVA guardada en: results/gsva_pathway_scores_by_sample.csv\n\n")

# 5. Evaluación Estadística basada en colData(vsd)
cond_vec <- as.character(metadata$condition)
ce_cols <- which(cond_vec %in% c("Control", "CE_MS", "CE"))
ee_cols <- which(cond_vec %in% c("Endometriosis", "EE_MS", "EE"))

if (length(ce_cols) == 0 || length(ee_cols) == 0) {
    ce_cols <- 1:(ncol(gsva_result)/2)
    ee_cols <- ((ncol(gsva_result)/2)+1):ncol(gsva_result)
}

stat_df <- data.frame(
    Pathway = rownames(gsva_result),
    Mean_Control = rowMeans(gsva_result[, ce_cols, drop = FALSE]),
    Mean_Endometriosis = rowMeans(gsva_result[, ee_cols, drop = FALSE]),
    p_value = apply(gsva_result, 1, function(x) {
        tryCatch(t.test(x[ee_cols], x[ce_cols])$p.value, error = function(e) NA)
    })
) %>%
    mutate(
        Delta_Score = Mean_Endometriosis - Mean_Control,
        Estado = ifelse(!is.na(p_value) & p_value < 0.05 & Delta_Score < 0, "Deprimida en Endometriosis",
                 ifelse(!is.na(p_value) & p_value < 0.05 & Delta_Score > 0, "Activada en Endometriosis", "Sin cambio sig."))
    )

write.csv(stat_df, "results/gsva_pathway_statistics.csv", row.names = FALSE)
cat("✓ Estadísticas comparativas guardadas en: results/gsva_pathway_statistics.csv\n\n")

# 6. Generación del Heatmap Anotado por Muestra
cat("Generando Heatmap anotado por muestra...\n")
annot_col <- data.frame(Condicion = cond_vec)
rownames(annot_col) <- colnames(gsva_result)

pheatmap(
    gsva_result,
    annotation_col = annot_col,
    main = "Actividad de Vías DDR por Muestra Individual (GSVA)",
    color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 8,
    filename = "plots/01_gsva_heatmap_annotated.png",
    width = 9,
    height = 6
)
cat("✓ Heatmap guardado en: plots/01_gsva_heatmap_annotated.png\n")

# 7. Generación del Boxplot para Reorganización Homóloga (HR)
cat("Generando Boxplot de Reorganización Homóloga (HR)...\n")
cond_labels <- ifelse(1:ncol(gsva_result) %in% ce_cols, "Control", "Endometriosis")

df_hr <- data.frame(
    Muestra = colnames(gsva_result),
    Condicion = factor(cond_labels, levels = c("Control", "Endometriosis")),
    Score = gsva_result["Homologous_Recombination", ]
)

p_box <- ggplot(df_hr, aes(x = Condicion, y = Score, fill = Condicion)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", width = 0.4) +
    geom_jitter(width = 0.12, size = 3.5, aes(color = Condicion)) +
    scale_fill_manual(values = c("Control" = "#2b5c8f", "Endometriosis" = "#d95f02")) +
    scale_color_manual(values = c("Control" = "#1a3856", "Endometriosis" = "#803801")) +
    theme_bw(base_size = 12) +
    labs(
        title = "Actividad de Reorganización Homóloga (HR)",
        subtitle = "Puntaje GSVA por Muestra Individual",
        x = "Grupo Biológico",
        y = "GSVA Score"
    ) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        legend.position = "none"
    )

ggsave("plots/02_gsva_hr_boxplot.png", p_box, width = 6, height = 5, dpi = 300)
cat("✓ Boxplot guardado en: plots/02_gsva_hr_boxplot.png\n\n")

cat("=================================================================\n")
cat("=== ANALISIS GSVA FINALIZADO EXITOSAMENTE =======================\n")
cat("=================================================================\n")
