#!/bin/bash
#SBATCH --job-name=qc_raw
#SBATCH --partition=internos
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/qc_raw_%j.out
#SBATCH --error=logs/qc_raw_%j.err

# Activar el entorno Conda
source ~/software/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate endometriosis_tpfinal

# Definir directorios
RAW_DIR="/home/ngs/Bioinfo3/TP_final/MondinoGasco/01_raw_data"
OUT_DIR="$(pwd)"

echo "=== Inicio de Control de Calidad (FastQC): $(date) ==="

# 1. Ejecutar FastQC usando --outdir
fastqc -t 4 --outdir ${OUT_DIR} ${RAW_DIR}/*.fastq.gz

# 2. Generar reporte integrado con MultiQC
echo "=== Generando reporte integrado con MultiQC ==="
multiqc ${OUT_DIR} -o ${OUT_DIR} -n multiqc_raw_report.html

echo "=== Control de Calidad finalizado: $(date) ==="
