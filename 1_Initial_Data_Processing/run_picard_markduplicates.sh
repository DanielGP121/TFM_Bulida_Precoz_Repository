#!/usr/bin/env bash
set -euo pipefail # Exit on error, undefined variable, or pipe failure
trap 'echo "ERROR: Script failed on line $LINENO with exit code $?"; exit 1' ERR # Error trapping

# run_picard_markduplicates.sh
# Marks PCR duplicates in coordinate-sorted BAM files using Picard MarkDuplicates.

# ── CONFIGURATION ────────────────────────────────────────────────────────────────
ALIGNDIR="/data/training2/analisis_TFM_Bulida_Precoz/02_alignment"
DEDUPDIR="/data/training2/analisis_TFM_Bulida_Precoz/03_deduplication_picard"
TMP_PICARD_DIR="${DEDUPDIR}/tmp_picard"
LOGDIR="${DEDUPDIR}/logs"
MEM_GB=20 # Memory in GB for Java Virtual Machine (Picard)

# ── PREPARATION ───────────────────────────────────────────────────────────────────
echo "INFO: Picard MarkDuplicates script started at $(date)"

if ! command -v picard >/dev/null; then
  echo "ERROR: picard command not found. Please ensure Picard Tools is installed and in your PATH (activate Conda env if needed)."
  exit 1
fi

if [[ ! -d "$ALIGNDIR" ]]; then
  echo "ERROR: Input alignment directory '$ALIGNDIR' does not exist. Run Stage 1 script first."
  exit 1
fi

mkdir -p "$DEDUPDIR" "$TMP_PICARD_DIR" "$LOGDIR"
if [[ ! -w "$DEDUPDIR" || ! -w "$TMP_PICARD_DIR" || ! -w "$LOGDIR" ]]; then
  echo "ERROR: Cannot write to one or more output/temp/log directories for deduplication."
  exit 1
fi

echo "===== SCRIPT CONFIGURATION (run_picard_markduplicates.sh) ====="
echo "Aligned BAM Dir         : $ALIGNDIR"
echo "Deduplicated Output Dir : $DEDUPDIR"
echo "Picard Temp Dir         : $TMP_PICARD_DIR"
echo "Logs Dir                : $LOGDIR"
echo "Java Memory (GB)        : $MEM_GB"
echo "============================================================"

#### MARK DUPLICATES LOOP ####
for SAMPLE in wt1 wt2 wt3 wt4 mt1 mt2 mt3 mt4; do
  echo "--- Marking duplicates for sample: $SAMPLE ---"

  INPUT_BAM="$ALIGNDIR/${SAMPLE}.aligned.sorted.bam"
  OUTPUT_DEDUP_BAM="$DEDUPDIR/${SAMPLE}.dedup.bam"
  METRICS_FILE="$DEDUPDIR/${SAMPLE}.dedup_metrics.txt"

  if [[ ! -f "$INPUT_BAM" ]]; then
    echo "WARNING: Input BAM file '$INPUT_BAM' not found for sample $SAMPLE. Skipping."
    continue
  fi
  if [[ ! -f "${INPUT_BAM}.bai" && ! -f "${INPUT_BAM%.bam}.bai" ]]; then
    echo "WARNING: Index file for '$INPUT_BAM' not found. Please ensure it is indexed. Skipping $SAMPLE."
    continue
  fi

  echo "Input BAM         : $INPUT_BAM"
  echo "Output Marked BAM : $OUTPUT_DEDUP_BAM"
  echo "Metrics File      : $METRICS_FILE"

  echo "Setting Java options and launching Picard MarkDuplicates for $SAMPLE..."

  # Set Java options using _JAVA_OPTIONS environment variable
  # This variable is often picked up by the JVM when launched via wrapper scripts.
  export _JAVA_OPTIONS="-Xmx${MEM_GB}G -Djava.io.tmpdir=${TMP_PICARD_DIR}"

  # Run Picard MarkDuplicates
  ( \
    picard MarkDuplicates \
    I="$INPUT_BAM" \
    O="$OUTPUT_DEDUP_BAM" \
    M="$METRICS_FILE" \
    REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR="$TMP_PICARD_DIR" \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    ASSUME_SORT_ORDER=coordinate \
    CREATE_INDEX=true && \
  echo "SUCCESS: Picard MarkDuplicates completed for $SAMPLE." \
  ) || \
  { echo "ERROR: Picard MarkDuplicates failed for $SAMPLE. Check script output/error logs."; }

  # Unset _JAVA_OPTIONS if you want to ensure it doesn't affect other Java applications
  unset _JAVA_OPTIONS
  echo "--- Finished marking duplicates for sample: $SAMPLE ---"
done

echo "============================================================"
echo "Duplicate marking finished for all applicable samples."
echo "BAM files with duplicates marked are in: $DEDUPDIR"
echo "============================================================"
echo "INFO: Script finished at $(date)"
