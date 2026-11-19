#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
})

cat("=== LÍNEA B - PASO 2: PRUEBA T DE STUDENT NO APAREADA EN FPKMS ===\n")

input_file <- "results/isoforms_fpkm_repair_genes.csv"
if (!file.exists(input_file)) {
    stop("ERROR: No se encontró '", input_file, "'.")
}

df <- read.csv(input_file)
genes <- unique(df$gene_symbol)
res_list <- list()

for (g in genes) {
    sub_df <- df %>% filter(gene_symbol == g)
    
    vals_ce <- sub_df %>% filter(condition == "Control") %>% pull(FPKM)
    vals_ee <- sub_df %>% filter(condition == "Endometriosis") %>% pull(FPKM)
    
    mean_ce <- mean(vals_ce, na.rm = TRUE)
    mean_ee <- mean(vals_ee, na.rm = TRUE)
    
    fc <- if (!is.na(mean_ce) && mean_ce > 0) mean_ee / mean_ce else NA
    log2fc <- if (!is.na(fc) && fc > 0) log2(fc) else NA
    
    p_val <- tryCatch({
        if (length(vals_ce) >= 2 && length(vals_ee) >= 2 && (sd(vals_ce) + sd(vals_ee)) > 0) {
            t.test(vals_ee, vals_ce, var.equal = FALSE)$p.value
        } else {
            NA
        }
    }, error = function(e) NA)
    
    res_list[[g]] <- data.frame(
        symbol = g,
        mean_FPKM_CE = mean_ce,
        mean_FPKM_EE = mean_ee,
        fold_change = fc,
        log2FoldChange = log2fc,
        p_value = p_val
    )
}

results_df <- bind_rows(res_list) %>% arrange(p_value)

write.csv(results_df, "results/t_test_fpkm_results.csv", row.names = FALSE)
cat("✓ Resultados de la prueba t de Student guardados en 'results/t_test_fpkm_results.csv'\n")
