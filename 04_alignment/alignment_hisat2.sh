#!/bin/bash
#SBATCH --job-name=hisat2_align
#SBATCH --partition=internos
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/align_%j.out
#SBATCH --error=logs/align_%j.err

# Activar el entorno Conda
source ~/software/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate endometriosis_tpfinal

# Directorios
BASE_DIR="/home/ngs/Bioinfo3/TP_final/MondinoGasco"
TRIM_DIR="${BASE_DIR}/03_trimming"
REF_INDEX="${BASE_DIR}/00_ref/grch38/genome"
OUT_DIR="$(pwd)"

echo "=== Inicio de Alineamiento con HISAT2: $(date) ==="

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
    echo "Alineando muestra: ${SRA}..."
    
    # HISAT2 -> SAM -> BAM ordenado
    hisat2 -p 8 \
        -x ${REF_INDEX} \
        -1 ${TRIM_DIR}/${SRA}_1.clean.fastq.gz \
        -2 ${TRIM_DIR}/${SRA}_2.clean.fastq.gz \
        --summary-file ${OUT_DIR}/${SRA}_hisat2_summary.txt \
        | samtools view -@ 4 -bS - \
        | samtools sort -@ 4 -o ${OUT_DIR}/${SRA}.sorted.bam

    # Indexar el archivo BAM
    samtools index ${OUT_DIR}/${SRA}.sorted.bam
done

# Generar informe consolidado de alineamiento con MultiQC
echo "=== Generando reporte MultiQC para Alineamiento ==="
multiqc ${OUT_DIR} -o ${OUT_DIR} -n multiqc_alignment_report.html

echo "=== Alineamiento finalizado: $(date) ==="
