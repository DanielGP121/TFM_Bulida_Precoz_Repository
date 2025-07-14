#!/usr/bin/env bash

# concat_mosdepth_summary.sh
# Script to concatenate mosdepth summary reports

# Directory where the mosdepth summary reports are located
REPORTS_DIR="/data/training2/analisis_TFM_Bulida_Precoz/02_alignment/coverage_stats"
# Consolidated output file
OUTPUT_FILE="${REPORTS_DIR}/all_samples.mosdepth_summary_combined.txt"

echo "Concatenating *.mosdepth.summary.txt files into ${OUTPUT_FILE}..."

# Iterate over all .mosdepth.summary.txt files in the specified directory
for f in "${REPORTS_DIR}"/*.mosdepth.summary.txt; do
  # Check if the file exists and is a regular file
  if [[ -f "$f" ]]; then
    echo "Processing: $f"
    # Add the filename as a header to the output file
    echo "=== Content of: $(basename "$f") ===" >> "$OUTPUT_FILE"
    # Add the content of the current file to the output file
    cat "$f" >> "$OUTPUT_FILE"
    # Add a couple of blank lines for better separation between file contents
    echo -e "\n\n" >> "$OUTPUT_FILE"
  fi
done

echo "Concatenation complete. The combined summary is in: $OUTPUT_FILE"
