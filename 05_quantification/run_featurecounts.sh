#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --partition=internos
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/quant_%j.out
#SBATCH --error=logs/quant_%j.err

# Activar entorno Conda
source ~/software/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate endometriosis_tpfinal

# Directorios
BASE_DIR="/home/ngs/Bioinfo3/TP_final/MondinoGasco"
ALIGN_DIR="${BASE_DIR}/04_alignment"
GTF_FILE="${BASE_DIR}/00_ref/Homo_sapiens.GRCh38.110.gtf"
OUT_DIR="$(pwd)"

echo "=== Inicio de Cuantificación con featureCounts: $(date) ==="

# 1. Ejecutar featureCounts para lecturas Paired-End
featureCounts -T 8 \
    -p -B -C \
    -t exon \
    -g gene_id \
    -a ${GTF_FILE} \
    -o ${OUT_DIR}/gene_counts_raw.txt \
    ${ALIGN_DIR}/*.sorted.bam

echo "=== Formateando matriz para R / DESeq2 ==="
# 2. Generar la matriz limpia de conteos (Geneid + 8 muestras)
# - Omitir comentario inicial (#)
# - Extraer columnas 1 y 7-14
# - Limpiar rutas y extensiones de los nombres de las muestras
grep -v "^#" ${OUT_DIR}/gene_counts_raw.txt \
    | cut -f1,7-14 \
    | sed 's|'${ALIGN_DIR}'/||g' \
    | sed 's/\.sorted\.bam//g' \
    > ${OUT_DIR}/gene_counts_matrix.txt

# 3. Generar la matriz de longitudes de genes (Geneid + Length)
grep -v "^#" ${OUT_DIR}/gene_counts_raw.txt \
    | cut -f1,6 \
    > ${OUT_DIR}/gene_lengths.txt

echo "=== Muestra de la matriz formateada ==="
head -n 5 ${OUT_DIR}/gene_counts_matrix.txt

echo "=== Generando reporte MultiQC de cuantificación ==="
multiqc ${OUT_DIR} -o ${OUT_DIR} -n multiqc_quant_report.html

echo "=== Cuantificación finalizada: $(date) ==="
