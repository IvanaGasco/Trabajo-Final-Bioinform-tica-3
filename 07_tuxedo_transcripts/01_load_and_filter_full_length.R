#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
})

cat("=== LÍNEA B - PASO 1: CARGA Y FILTRADO DE TRANSCRITOS FULL-LENGTH ===\n")

repair_genes <- list(
  MMR = c("MLH1", "MSH3", "MSH6", "PCNA", "EXO1", "LIG1"),
  BER = c("OGG1", "UNG", "LIG1", "FEN1", "POLB", "PNKP", "RFC1", "PCNA", "SMUG1", "MBD4", "TDG", "NEIL3", "PARP2", "APEX1"),
  NER = c("DDB2", "RFC1", "PCNA", "GTF2H6", "LIG1", "FEN1", "CETN2", "RPA2", "RPA3", "GTF2H2", "MNAT1"),
  HR = c("RAD51", "DMC1", "RAD54B", "RAD54L", "BLM", "BRCA2", "BARD1", "MRE11", "EXO1", "GEN1", "SLX1B", "RAD51B", "BRCA1", "EME2", "ATM"),
  FA = c("FANCD2", "FANCI", "BRIP1", "FANCO", "FAAP20", "FAAP24", "BRCA2"),
  NHEJ = c("XRCC4", "DCLRE1C", "PNKP", "POLM", "POLL")
)

target_genes <- unique(unlist(repair_genes))

fpkm_long_path <- "results/isoforms_fpkm_long.csv"
if (!file.exists(fpkm_long_path)) {
    stop("ERROR: No se encontró '", fpkm_long_path, "'.")
}

fpkm_long <- read.csv(fpkm_long_path)

fpkm_filtered <- fpkm_long %>%
    filter(gene_symbol %in% target_genes) %>%
    mutate(condition = ifelse(sample %in% c("SRR12137399", "SRR12137400", "SRR12137401", "SRR12137402"), "Control", "Endometriosis"))

write.csv(fpkm_filtered, "results/isoforms_fpkm_repair_genes.csv", row.names = FALSE)
cat("✓ Transcritos filtrados para genes de reparación guardados en 'results/isoforms_fpkm_repair_genes.csv'\n")
