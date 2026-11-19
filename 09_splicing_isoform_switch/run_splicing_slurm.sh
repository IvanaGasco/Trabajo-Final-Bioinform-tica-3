#!/bin/bash
#SBATCH --job-name=isoform_switch
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

echo "=== Ejecutando Análisis de Splicing en SLURM ==="
echo "Nodo: $(hostname)"
echo "Fecha: $(date)"

# Crear directorio de logs si no existe
mkdir -p logs

# Cargar entorno de R/Conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate endometriosis_tpfinal

Rscript 01_isoform_switch_analysis.R

echo "=== Finalizado: $(date) ==="
