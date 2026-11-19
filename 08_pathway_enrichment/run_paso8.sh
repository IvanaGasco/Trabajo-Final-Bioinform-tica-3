#!/bin/bash
#SBATCH --job-name=enrichment_paso8
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

# Activar el entorno Conda
eval "$(conda shell.bash hook)"
conda activate endometriosis_tpfinal

echo "=== EJECUTANDO PASO 8: ENRIQUECIMIENTO FUNCIONAL EN SLURM ==="
Rscript 01_pathway_enrichment.R
