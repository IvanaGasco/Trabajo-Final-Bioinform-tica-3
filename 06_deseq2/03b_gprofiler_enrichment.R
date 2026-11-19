#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(gprofiler2)
    library(dplyr)
    library(ggplot2)
})

cat("=== SCRIPT 3B: ENRIQUECIMIENTO FUNCIONAL DINÁMICO (g:Profiler) ===\n")

dir.create("results_gprofiler", showWarnings = FALSE, recursive = TRUE)

degs_path <- "results/degs_significant_FDR0.05_FC1.5.csv"
if (!file.exists(degs_path)) {
    stop("ERROR: No se encontró '", degs_path, "'. Ejecutá primero el pipeline principal en SLURM.")
}

degs <- read.csv(degs_path)
if (nrow(degs) == 0) {
    stop("ERROR: La lista de genes significativos está vacía.")
}

ensembl_ids <- sub("\\..*", "", degs$gene_id)

cat("Consultando API web de g:Profiler para", length(ensembl_ids), "genes Ensembl...\n")

gostres <- gost(
    query = ensembl_ids,
    organism = "hsapiens",
    sources = c("GO:BP", "KEGG", "REAC")
)

if (!is.null(gostres$result)) {
    res_clean <- gostres$result %>%
        select(source, term_id, term_name, p_value, term_size, query_size, intersection_size) %>%
        arrange(p_value)

    out_csv <- "results_gprofiler/gprofiler2_enrichment_results.csv"
    write.csv(res_clean, out_csv, row.names = FALSE)
    cat("✓ Tabla guardada dinámicamente en:", out_csv, "(", nrow(res_clean), "vías encontradas)\n")

    out_png <- "results_gprofiler/06_gprofiler2_enrichment_plot.png"
    p_gost <- gostplot(gostres, capped = TRUE, interactive = FALSE)
    ggsave(out_png, p_gost, width = 10, height = 6, dpi = 300)
    cat("✓ Gráfico guardado en:", out_png, "\n")
} else {
    cat("Atención: No se encontraron términos enriquecidos significativos.\n")
}

cat("=== SCRIPT 3B FINALIZADO CON ÉXITO ===\n")
