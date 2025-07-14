#!/usr/bin/env bash

#concat_flagstats.sh
# Script to concatenate flagstat reports
# Directory where the reports are located
REPORTS_DIR="/data/training2/analisis_TFM_Bulida_Precoz/02_alignment/flagstat_reports"
# Consolidated output file
OUTPUT_FILE="${REPORTS_DIR}/all_samples.flagstat_summary.txt"

echo "Concatenating .flagstat.txt files into ${OUTPUT_FILE}..."
# Iterate over all .flagstat.txt files in the specified directory
for f in "${REPORTS_DIR}"/*.flagstat.txt; do
  # Check if the file exists (the glob might not find anything)
  if [[ -f "$f" ]]; then
    echo "Processing: $f"
    # Add the filename as a header to the output file
    echo "=== Content of: $(basename "$f") ===" >> "$OUTPUT_FILE"
    # Add the content of the current file to the output file
    cat "$f" >> "$OUTPUT_FILE"
    # Add a blank line (or two) to separate sections
    echo -e "\n" >> "$OUTPUT_FILE"
  fi
done

echo "Concatenation complete. The summary is in: $OUTPUT_FILE"
