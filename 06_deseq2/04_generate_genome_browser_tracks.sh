#!/bin/bash
# ==============================================================================
# SCRIPT 4: GENERACIÓN DE PISTAS PARA GENOME BROWSER
# ==============================================================================

echo "=== SCRIPT 4: GENERACIÓN DE PISTAS PARA GENOME BROWSER ==="

mkdir -p results_genome_browser

# Archivo BED con coordenadas de interés para los genes clave del trabajo
cat << 'BEDEOF' > results_genome_browser/genes_interes_igv.bed
chr24778300047807000MSH6
chr154069500040732000RAD51
chr165662000056625000MT1G
chr98701000087030000GADD45G
BEDEOF

echo "✓ Archivo BED de genes clave creado en results_genome_browser/genes_interes_igv.bed"

BAM_DIR=""
for d in "../05_alignment" "../04_alignment" "alignment"; do
    if [ -d "$d" ]; then
        BAM_DIR="$d"
        break
    fi
done

if [ -z "$BAM_DIR" ]; then
    BAM_DIR="../05_alignment"
fi

found_bams=0
for bam in "$BAM_DIR"/*.sorted.bam "$BAM_DIR"/*.bam; do
    if [ -f "$bam" ]; then
        found_bams=1
        sample=$(basename "$bam" | sed 's/\.sorted\.bam//' | sed 's/\.bam//')
        echo "Procesando $sample..."

        # Indexar si falta el archivo .bai
        if [ ! -f "${bam}.bai" ] && [ ! -f "${BAM_DIR}/${sample}.sorted.bam.bai" ] && [ ! -f "${BAM_DIR}/${sample}.bam.bai" ]; then
            samtools index "$bam"
        fi

        # Generar pista BigWig con bamCoverage si está disponible
        if command -v bamCoverage &> /dev/null; then
            bamCoverage -b "$bam" -o "results_genome_browser/${sample}.bw" --normalizeUsing CPM -p 4
        else
            echo "Aviso: bamCoverage no disponible, copiando BAM/BAI..."
            cp "$bam" "results_genome_browser/${sample}.bam"
            if [ -f "${bam}.bai" ]; then
                cp "${bam}.bai" "results_genome_browser/${sample}.bam.bai"
            elif [ -f "${bam}.sorted.bam.bai" ]; then
                cp "${bam}.sorted.bam.bai" "results_genome_browser/${sample}.bam.bai"
            fi
        fi
    fi
done

if [ $found_bams -eq 0 ]; then
    echo "Atención: No se encontraron archivos BAM en $BAM_DIR."
else
    echo "Script 4 finalizado con éxito."
fi
