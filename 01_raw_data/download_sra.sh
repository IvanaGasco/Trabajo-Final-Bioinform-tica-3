#!/bin/bash
#SBATCH --job-name=ena_download
#SBATCH --partition=internos
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err

OUT_DIR="$(pwd)"
echo "=== Inicio de descarga desde ENA: $(date) ==="

# wget -c permite reanudar descargas si se interrumpe la red
wget -c -i urls.txt -P ${OUT_DIR}

echo "=== Descarga finalizada: $(date) ==="
