# Trabajo Práctico Final - Bioinformática 3 (FIUNER)##  Proyecto: Modulación de la Respuesta al Daño en el ADN (DDR) en Endometrio
- **Autoras:** Mondino & Gascó
- **Materia:** Bioinformática 3 - Facultad de Ingeniería (FIUNER)
- **Dataset Principal:** GSE153740 (*Homo sapiens*, Bulk RNA-seq Paired-End 2x76 bp, dUTP strand-specific)
- **Diseño Experimental:** 8 muestras en Fase Secretora Media (4 Controles `CE-MS` vs. 4 Endometriosis `EE-MS`)
- **Entorno Conda:** `endometriosis_tpfinal` (Python 3.10 / R 4.x)

---## Creación y Configuración del Entorno Conda

# 1. Creación del entorno
conda create -n endometriosis_tpfinal python=3.10 -y
# 2. Configuración estricta de canales
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

## Lista de Herramientas y Paquetes Instalados 
### 1. Herramientas de Línea de Comandos (CLI / Bioconda & Conda-Forge)
- **sra-tools** (bioconda):
  - **Uso:** [Paso 1] Obtención de datos crudos desde NCBI SRA mediante prefetch y extracción de archivos FASTQ pareados con fasterq-dump para las 8 muestras del estudio (SRR12137399 a SRR12137406).
- **fastqc** (bioconda):
  - **Uso:** [Paso 2] Control de calidad inicial (lecturas crudas) y final (post-trimming), evaluando perfiles de calidad Phred (Q), contenido GC, adaptadores e índices de duplicación.
- **fastp** (bioconda):
  - **Uso:** [Paso 3] Preprocesamiento ultra-rápido, recorte de adaptadores de biblioteca Illumina, eliminación de bases con baja calidad (Q < 20) y filtrado de lecturas por longitud mínima.
- **multiqc** (bioconda / conda-forge):
  - **Uso:** [Paso 2] Integración y consolidación de los informes individuales de FastQC, fastp, HISAT2 y featureCounts en un único reporte HTML interactivo.
- **hisat2** (bioconda):
  - **Uso:** [Paso 4] Alineamiento genómico splice-aware de lecturas limpias contra el genoma de referencia Homo sapiens (GRCh38), configurado con orientación específica de hebra dUTP (`--rna-strandness RF`).
- **samtools** (bioconda):
  - **Uso:** [Paso 4] Procesamiento de alineamientos, conversión SAM a BAM, ordenamiento por coordenadas genómicas (`samtools sort`) e indexación (`samtools index`).
- **subread / featureCounts** (bioconda):
  - **Uso:** [Paso 5] Cuantificación de expresión a nivel de gen basada en la anotación GTF de Ensembl, generando la matriz de recuentos crudos para la Línea A (Pipeline Cátedra).
- **stringtie** (bioconda):
  - **Uso:** [Paso 7] Cuantificación y ensamblado de abundancia a nivel de isoformas/transcritos en FPKM/TPM para la Línea B (Pipeline Paper Original) y entrada al análisis de splicing.
- **deeptools** (bioconda):
  - **Uso:** [Paso 4 / Paso 10] Generación y normalización de pistas de cobertura genómica BigWig (`.bw`) mediante `bamCoverage` para inspección visual de densidad de lecturas exónicas en IGV Browser.

### 2. Paquetes de R y Bioconductor (Análisis Estadístico, Visualización y Minería)
- **r-base**: Entorno de programación estadística.
- **bioconductor-deseq2** (DESeq2):
  - **Uso:** [Paso 6] Normalización de conteos por gen, transformación de estabilización de varianza (VST), análisis de componentes principales (PCA), matriz de distancias y análisis de expresión diferencial (DEGs).
- **bioconductor-gsva** (GSVA):
  - **Uso:** [Paso 10] Cálculo de puntajes de actividad por muestra individual (Gene Set Variation Analysis / ssGSEA) sobre firmas biológicas curadas de la vía DDR (HR, BER, NER, MMR, NHEJ, FA, p53).
- **r-gprofiler2** (gprofiler2):
  - **Uso:** [Paso 8] Análisis de enriquecimiento funcional comparativo en ontologías y vías (`GO:BP`, `KEGG`, `Reactome`, `WikiPathways`) entre los DEGs de la Línea A y los transcritos alterados de la Línea B.
- **r-pheatmap** (pheatmap):
  - **Uso:** [Paso 6 / Paso 10] Construcción de Heatmaps con agrupamiento jerárquico (clustering) y anotación clínica de condiciones (Controles vs. Endometriosis).
- **r-ggplot2** (ggplot2):
  - **Uso:** [Paso 6 / Puntos 8, 9, 10] Generación de gráficos científicos de alta resolución: Volcano plots, MA plots, Barplots de proporciones de isoformas y Boxplots comparativos de actividad por vía.
- **r-tidyverse** (dplyr, tidyr, readr):
  - **Uso:** [Pasos 7, 8, 9, 10] Limpieza, filtrado, reestructuración de dataframes de expresión e isoformas, y cálculo de la Fracción de Isoforma ($\text{IF}$) y DeltaIF.
- **r-rcolorbrewer** (RColorBrewer):
  - **Uso:** [Pasos 6, 8, 9, 10] Paletas de colores profesionales y accesibles para heatmaps y figuras del manuscrito.
