#!/bin/bash
#SBATCH --job-name=quant_transcripts
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

set -e

cd /home/ngs/Bioinfo3/TP_final/MondinoGasco/07_tuxedo_transcripts

source ~/.bashrc 2>/dev/null || true
conda activate endometriosis_tpfinal 2>/dev/null || true

mkdir -p logs results FPKM_per_sample

echo "=== LÍNEA B: CUANTIFICACIÓN DE FPKM A NIVEL DE TRANSCRITO ==="

# 1. Buscar archivo GTF de referencia
GTF_PATH=""
for p in ../00_annotation/*.gtf ../00_genome/*.gtf ../05_quantification/*.gtf; do
    if [ -f "$p" ]; then
        GTF_PATH="$p"
        break
    fi
done

if [ -z "$GTF_PATH" ]; then
    GTF_PATH=$(find /home/ngs/Bioinfo3/TP_final/MondinoGasco/ -name "*.gtf" | head -n 1)
fi

echo "Anotación GTF detectada: $GTF_PATH"

# 2. Procesar las 8 muestras
SAMPLES=("SRR12137399" "SRR12137400" "SRR12137401" "SRR12137402" "SRR12137403" "SRR12137404" "SRR12137405" "SRR12137406")

for s in "${SAMPLES[@]}"; do
    BAM="../04_alignment/${s}.sorted.bam"
    OUT_DIR="FPKM_per_sample/${s}"
    mkdir -p "$OUT_DIR"
    
    echo "Procesando $s a nivel transcrito..."
    if command -v stringtie &> /dev/null; then
        stringtie "$BAM" -G "$GTF_PATH" -e -B -p 4 -A "$OUT_DIR/gene_abund.tab" -o "$OUT_DIR/transcripts.gtf"
    elif command -v cufflinks &> /dev/null; then
        cufflinks -p 4 -G "$GTF_PATH" -o "$OUT_DIR" "$BAM"
    fi
done

# 3. Consolidar matriz en R
Rscript - << 'RSCRIPT'
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

samples <- c("SRR12137399", "SRR12137400", "SRR12137401", "SRR12137402", "SRR12137403", "SRR12137404", "SRR12137405", "SRR12137406")
fpkm_list <- list()

for (s in samples) {
    # Buscar t_data.ctab o isoforms.fpkm_tracking o stringtie output
    ctab_path <- file.path("FPKM_per_sample", s, "t_data.ctab")
    cuff_path <- file.path("FPKM_per_sample", s, "isoforms.fpkm_tracking")
    
    if (file.exists(ctab_path)) {
        df <- read.table(ctab_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        sub_df <- df %>% select(t_name, gene_name, FPKM) %>% rename(transcript_id = t_name, gene_symbol = gene_name)
    } else if (file.exists(cuff_path)) {
        df <- read.table(cuff_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        sub_df <- df %>% select(tracking_id, gene_short_name, FPKM) %>% rename(transcript_id = tracking_id, gene_symbol = gene_short_name)
    } else {
        next
    }
    sub_df$sample <- s
    fpkm_list[[s]] <- sub_df
}

if (length(fpkm_list) > 0) {
    all_fpkm <- bind_rows(fpkm_list)
    write.csv(all_fpkm, "results/isoforms_fpkm_long.csv", row.names = FALSE)
    
    matrix_fpkm <- all_fpkm %>%
        pivot_wider(id_cols = c(transcript_id, gene_symbol), names_from = sample, values_from = FPKM, values_fill = 0)
        
    write.csv(matrix_fpkm, "results/isoforms_fpkm_matrix.csv", row.names = FALSE)
    cat("✓ Matriz de FPKMs consolidada en 'results/isoforms_fpkm_matrix.csv' (", nrow(matrix_fpkm), "transcritos)\n")
} else {
    cat("Advertencia: No se encontraron archivos de salida de cuantificación.\n")
}
RSCRIPT

echo "=== PASO 1 LÍNEA B FINALIZADO EXITOSAMENTE ==="
