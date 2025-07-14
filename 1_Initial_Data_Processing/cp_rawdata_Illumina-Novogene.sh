#!/usr/bin/env bash
# Exit immediately if a command exits with a non-zero status (-e),
# treat unset variables as an error (-u),
# and ensure any failure in a pipeline causes the whole pipeline to fail (-o pipefail)
set -euo pipefail

# cp_rawdata_Illumina-Novogene.sh
# Copy all Illumina-Novogene FASTQ files into the working directory
#
# RAWDIR: location of the raw_data subfolders
# OUTDIR: destination for all .fq.gz files

RAWDIR="/data/training2/reads/seqs/Illumina-Novogene/X204SC20112838-Z01-F001/raw_data"
OUTDIR="/data/training2/analisis_TFM_Bulida_Precoz/00_raw_data"

# Quick sanity‐check: print script name and key variables
echo "===== DEBUG: verifying script and paths ====="
echo "Script: $0"
echo "RAWDIR: $RAWDIR"
echo "OUTDIR : $OUTDIR"
echo "============================================="

echo "Starting copy of raw Illumina-Novogene FASTQ files..."

for fichero in "$RAWDIR"/*/*.fq.gz; do
    echo "Copying: cp $fichero $OUTDIR/"
    cp "$fichero" "$OUTDIR/"
done

echo "Copy of Illumina-Novogene FASTQ files completed."
