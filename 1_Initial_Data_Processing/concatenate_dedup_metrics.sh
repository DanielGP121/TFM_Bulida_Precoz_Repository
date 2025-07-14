#!/usr/bin/env bash

# concatenate_dedup_metrics.sh
# Script to concatenate Picard MarkDuplicates metrics files

# Directory where the MarkDuplicates metrics files are located
METRICS_DIR="/data/training2/analisis_TFM_Bulida_Precoz/03_deduplication_picard" 
# Consolidated output file
OUTPUT_FILE="${METRICS_DIR}/all_samples.dedup_metrics_combined.txt"

# Ensure the output file is empty at the start
> "$OUTPUT_FILE"

echo "Concatenating *.dedup_metrics.txt files into ${OUTPUT_FILE}..."

# Iterate over all .dedup_metrics.txt files in the specified directory
for f in "${METRICS_DIR}"/*.dedup_metrics.txt; do
  # Check if the file exists and is a regular file
  if [[ -f "$f" ]]; then
    echo "Processing: $f"
    # Add the filename as a header to the output file
    # $(basename "$f") extracts just the filename from the full path
    echo "=== Metrics from: $(basename "$f") ===" >> "$OUTPUT_FILE"
    # Add the content of the current file to the output file
    cat "$f" >> "$OUTPUT_FILE"
    # Add a couple of blank lines for better separation between file contents
    echo -e "\n\n" >> "$OUTPUT_FILE"
  fi
done

echo "Concatenation complete. The combined metrics are in: $OUTPUT_FILE"
