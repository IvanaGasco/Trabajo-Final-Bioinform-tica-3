#!/bin/bash

cd /home/ngs/Bioinfo3/TP_final/MondinoGasco/06_deseq2

source ~/.bashrc 2>/dev/null || true
conda activate endometriosis_tpfinal 2>/dev/null || true

Rscript 03b_gprofiler_enrichment.R
