#!/bin/bash
#SBATCH --job-name=deeptools_bw
#SBATCH --partition=internos
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/bigwig_%j.out
#SBATCH --error=logs/bigwig_%j.err

# Activar el entorno Conda
source ~/software/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate endometriosis_tpfinal

# Directorio de trabajo
ALIGN_DIR="$(pwd)"

echo "=== Inicio de generación de BigWig con bamCoverage: $(date) ==="

SAMPLES=(
    "SRR12137399"
    "SRR12137400"
    "SRR12137401"
    "SRR12137402"
    "SRR12137403"
    "SRR12137404"
    "SRR12137405"
    "SRR12137406"
)

for SRA in "${SAMPLES[@]}"; do
    echo "Generando BigWig para: ${SRA}..."
    bamCoverage \
        -b ${ALIGN_DIR}/${SRA}.sorted.bam \
        -o ${ALIGN_DIR}/${SRA}.bw \
        --normalizeUsing CPM \
        --binSize 10 \
        -p 8
done

echo "=== Generación de BigWig finalizada: $(date) ==="
