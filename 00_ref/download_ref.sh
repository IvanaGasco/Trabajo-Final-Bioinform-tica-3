#!/bin/bash
#SBATCH --job-name=get_ref
#SBATCH --partition=internos
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err

OUT_DIR="$(pwd)"
echo "=== Inicio de descarga: $(date) ==="

# 1. Descargar e descomprimir genoma GRCh38
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz -O ${OUT_DIR}/grch38_genome.tar.gz
tar -xzvf ${OUT_DIR}/grch38_genome.tar.gz -C ${OUT_DIR}
rm ${OUT_DIR}/grch38_genome.tar.gz

# 2. Descargar e descomprimir anotación GTF
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz -O ${OUT_DIR}/Homo_sapiens.GRCh38.110.gtf.gz
gunzip ${OUT_DIR}/Homo_sapiens.GRCh38.110.gtf.gz

echo "=== Descargas finalizadas: $(date) ==="
