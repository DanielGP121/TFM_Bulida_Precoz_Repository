#!/usr/bin/env bash
# SCRIPT: Runs FastQC on all raw read files to generate quality control reports.

set -euo pipefail

## CONFIGURATION ##
BASE_DIR="/data/training2/analisis_TFM_Bulida_Precoz"
# Directory containing the raw sequencing data
RAW_DATA_DIR="${BASE_DIR}/00_raw_data"
# A new directory to store the FastQC reports for this analysis
OUTPUT_DIR="${BASE_DIR}/09_samples_double_check/fastqc_reports"
mkdir -p "${OUTPUT_DIR}"

# Number of threads to use for FastQC
THREADS=8

## SCRIPT LOGIC ##
echo "--- Starting FastQC Analysis for All Raw Samples ---"

# Check if the raw data directory exists
if [ ! -d "${RAW_DATA_DIR}" ]; then
    echo "ERROR: Raw data directory not found at ${RAW_DATA_DIR}"
    exit 1
fi

# Find all fastq files (both .fastq.gz and .fq.gz) and run FastQC on them
# The output reports will be placed in the specified OUTPUT_DIR
fastqc \
    --threads "${THREADS}" \
    --outdir "${OUTPUT_DIR}" \
    ${RAW_DATA_DIR}/*.fastq.gz \
    ${RAW_DATA_DIR}/*.fq.gz

echo -e "\n### FastQC analysis complete! ###"
echo "Reports have been generated in: ${OUTPUT_DIR}"