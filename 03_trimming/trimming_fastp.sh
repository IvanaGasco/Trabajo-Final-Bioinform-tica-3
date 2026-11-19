#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --partition=internos
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/trim_%j.out
#SBATCH --error=logs/trim_%j.err

# Activar el entorno Conda
source ~/software/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate endometriosis_tpfinal

# Directorios
RAW_DIR="/home/ngs/Bioinfo3/TP_final/MondinoGasco/01_raw_data"
OUT_DIR="$(pwd)"

echo "=== Inicio de Trimming con fastp: $(date) ==="

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
    echo "Procesando muestra: ${SRA}..."
    fastp \
        --in1 ${RAW_DIR}/${SRA}_1.fastq.gz \
        --in2 ${RAW_DIR}/${SRA}_2.fastq.gz \
        --out1 ${OUT_DIR}/${SRA}_1.clean.fastq.gz \
        --out2 ${OUT_DIR}/${SRA}_2.clean.fastq.gz \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 35 \
        --thread 4 \
        --html ${OUT_DIR}/${SRA}_fastp.html \
        --json ${OUT_DIR}/${SRA}_json.json
done

# Generar reporte MultiQC post-trimming
echo "=== Generando reporte MultiQC post-trimming ==="
multiqc ${OUT_DIR} -o ${OUT_DIR} -n multiqc_trim_report.html

echo "=== Trimming finalizado: $(date) ==="
