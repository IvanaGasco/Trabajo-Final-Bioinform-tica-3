#!/bin/bash
#SBATCH --job-name=line_b_tuxedo
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

set -e

cd /home/ngs/Bioinfo3/TP_final/MondinoGasco/07_tuxedo_transcripts

source ~/.bashrc 2>/dev/null || true
conda activate endometriosis_tpfinal 2>/dev/null || true

mkdir -p logs results plots

echo "=== EJECUTANDO PRUEBA T DE STUDENT Y VISUALIZACIONES (LÍNEA B) ==="

Rscript 01_load_and_filter_full_length.R
Rscript 02_student_t_test_fpkm.R
Rscript 03_lineB_visualizations.R

echo "=== LÍNEA B FINALIZADA EXITOSAMENTE ==="
