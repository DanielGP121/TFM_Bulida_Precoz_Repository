#!/usr/bin/env bash

set -euo pipefail  
# Exit on error, undefined var or failed pipe

# concat_reads.sh

# Concatenate Illumina reads for samples E–H from projects 3649 & 4024
# Usage:
#   nohup bash scripts/concat_reads.sh > concat_reads.log 2>&1 &

# RAWDIR: base directory with all read subfolders
# OUTDIR: where to write the combined FASTQ files

RAWDIR="/data/training2/reads/seqs"
OUTDIR="/data/training2/analisis_TFM_Bulida_Precoz/00_raw_data"

# Quick sanity‐check: print script name and key variables
echo "===== DEBUG: verify environment ====="
echo "Script:    $0"
echo "RAWDIR:    $RAWDIR"
echo "OUTDIR:    $OUTDIR"
echo "====================================="

echo "Starting concatenation of samples E–H..."

# Concatenation loop: zcat → gzip
for SAMPLE in E F G H; do
  for READ in R1 R2; do
    echo "Processing sample ${SAMPLE}, read ${READ}…"
    zcat \
      "$RAWDIR"/3649_"${SAMPLE}"_*_"${READ}"_001.fastq.gz \
      "$RAWDIR"/4024/run562/4024_"${SAMPLE}"_*_"${READ}"_001.fastq.gz \
      "$RAWDIR"/4024/run564/4024_"${SAMPLE}"_*_"${READ}"_001.fastq.gz \
    | gzip > "$OUTDIR"/"${SAMPLE}_${READ}_combined.fastq.gz"
  done
done

echo "All concatenations completed successfully"
