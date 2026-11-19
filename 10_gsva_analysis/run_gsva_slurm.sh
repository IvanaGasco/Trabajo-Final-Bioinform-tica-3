#!/bin/bash
#SBATCH --job-name=gsva_paso10
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

echo "=== Ejecutando Paso 10: GSVA & ssGSEA en SLURM ==="
echo "Nodo asignado: $SLURMD_NODENAME"
echo "Fecha de inicio: $(date)"
echo "--------------------------------------------------"

# Asegurar que exista el directorio de logs
mkdir -p logs

# Cargar entorno de conda de forma robusta
source $(conda info --base)/etc/profile.d/conda.sh
conda activate endometriosis_tpfinal

Rscript 01_gsva_ssgsea_analysis.R

echo "--------------------------------------------------"
echo "Fecha de finalización: $(date)"
echo "=== Trabajo Finalizado ==="
