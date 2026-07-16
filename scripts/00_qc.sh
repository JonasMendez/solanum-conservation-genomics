#!/usr/bin/env bash
# =============================================================================
# 00_qc.sh
# Run FastQC + MultiQC on all demultiplexed FASTQs.
# Conda env: qc  (fastqc, multiqc)
# =============================================================================
set -euo pipefail

FASTQ_DIR="fastq"
OUT="00_qc"
THREADS=8

mkdir -p "${OUT}/fastqc"

echo "Running FastQC on $(ls ${FASTQ_DIR}/*.fastq.gz | wc -l) files..."
fastqc -t "${THREADS}" -o "${OUT}/fastqc" "${FASTQ_DIR}"/*.fastq.gz

echo "Running MultiQC..."
multiqc -o "${OUT}" "${OUT}/fastqc"

echo ""
echo "QC complete. Review: ${OUT}/multiqc_report.html"
echo "Check for: adapter content, low base quality, anomalous duplication rates"
echo "Flag any samples that look like outliers before proceeding to assembly."
