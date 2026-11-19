#!/bin/bash
#SBATCH --job-name=deseq2_pipeline
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

set -e

cd /home/ngs/Bioinfo3/TP_final/MondinoGasco/06_deseq2

source ~/.bashrc 2>/dev/null || true
conda activate endometriosis_tpfinal 2>/dev/null || true

mkdir -p logs plots results results_gprofiler results_genome_browser

echo "=== INICIANDO PIPELINE EN SLURM ==="

Rscript 01_exploratory_analysis.R
Rscript 02_differential_expression.R
Rscript 03_visualizations_dna_repair.R
bash 04_generate_genome_browser_tracks.sh

echo "=== PIPELINE EN SLURM FINALIZADO EXITOSAMENTE ==="
