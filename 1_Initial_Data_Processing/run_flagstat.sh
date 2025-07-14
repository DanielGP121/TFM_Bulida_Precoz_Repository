#!/usr/bin/env bash
set -euo pipefail
trap 'echo "ERROR: Script failed on line $LINENO with exit code $?"; exit 1' ERR

# run_flagstat.sh
# Calculates alignment statistics using samtools flagstat for all processed BAM files.

# ── CONFIGURATION ────────────────────────────────────────────────────────────────
ALIGNDIR="/data/training2/analisis_TFM_Bulida_Precoz/02_alignment"
FLAGSTAT_DIR="${ALIGNDIR}/flagstat_reports" # Directory to store flagstat reports

# ── PREPARE ────────────────────────────────────────────────────────────────────
echo "INFO: Script started at $(date)"

if ! command -v samtools >/dev/null; then
  echo "ERROR: samtools command not found. Activate Conda environment or install."; exit 1
fi

if [[ ! -d "$ALIGNDIR" ]]; then
  echo "ERROR: Alignment directory '$ALIGNDIR' does not exist. Run alignment script first."; exit 1
fi

mkdir -p "$FLAGSTAT_DIR"
if [[ ! -w "$FLAGSTAT_DIR" ]]; then
  echo "ERROR: Cannot write to flagstat reports directory '$FLAGSTAT_DIR'"; exit 1
fi

echo "===== SCRIPT CONFIGURATION (run_flagstat.sh) ====="
echo "Aligned BAM Dir      : $ALIGNDIR"
echo "Flagstat Reports Dir : $FLAGSTAT_DIR"
echo "================================================="

#### FLAGSTAT LOOP ####
for SAMPLE in wt1 wt2 wt3 wt4 mt1 mt2 mt3 mt4; do
  echo "--- Calculating flagstat for sample $SAMPLE ---"

  INPUT_BAM="$ALIGNDIR/${SAMPLE}.aligned.sorted.bam"
  OUTPUT_FLAGSTAT_FILE="$FLAGSTAT_DIR/${SAMPLE}.flagstat.txt"

  # Check if input BAM file exists
  if [[ ! -f "$INPUT_BAM" ]]; then
    echo "WARNING: Aligned BAM file '$INPUT_BAM' not found for sample $SAMPLE. Skipping."
    continue
  fi

  echo "Input BAM: $INPUT_BAM"
  echo "Output File: $OUTPUT_FLAGSTAT_FILE"

  # Run samtools flagstat and redirect output to a file
  echo "Running samtools flagstat for $SAMPLE..."
  samtools flagstat "$INPUT_BAM" > "$OUTPUT_FLAGSTAT_FILE"
  if [ $? -eq 0 ]; then
    echo "SUCCESS: samtools flagstat completed for $SAMPLE. Report saved to $OUTPUT_FLAGSTAT_FILE"
  else
    echo "ERROR: samtools flagstat failed for $SAMPLE."
    continue
  fi

  echo "--- Finished flagstat for sample $SAMPLE ---"
done

echo "================================================="
echo "samtools flagstat processing finished for all samples."
echo "Reports are in: $FLAGSTAT_DIR"
echo "================================================="
echo "INFO: Script finished at $(date)"
