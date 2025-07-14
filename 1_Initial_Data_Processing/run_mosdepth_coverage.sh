#!/usr/bin/env bash
set -euo pipefail
trap 'echo "ERROR: Script failed on line $LINENO with exit code $?"; exit 1' ERR

# run_mosdepth_coverage.sh
# Calculates depth of coverage using mosdepth for all aligned BAM files.

# ── CONFIGURATION ────────────────────────────────────────────────────────────────
ALIGNDIR="/data/training2/analisis_TFM_Bulida_Precoz/02_alignment"
COVERAGE_BASE_DIR="/data/training2/analisis_TFM_Bulida_Precoz/02_alignment/coverage_stats" # Base directory for all coverage outputs
REF_GENOME_FASTA="/data/training2/info/assemblies/BUL_cur_guided.v1.0.fasta" # Reference FASTA file
THREADS=8 # Number of threads to use for mosdepth

# ── PREPARE ────────────────────────────────────────────────────────────────────
echo "INFO: Script started at $(date)"

if ! command -v mosdepth >/dev/null; then
  echo "ERROR: mosdepth command not found. Please ensure it is installed and your Conda environment is activated."
  exit 1
fi
if ! command -v samtools >/dev/null; then # samtools needed for faidx check
  echo "ERROR: samtools command not found. Please ensure it is installed."
  exit 1
fi

if [[ ! -d "$ALIGNDIR" ]]; then
  echo "ERROR: Alignment directory '$ALIGNDIR' does not exist. Run alignment script first."
  exit 1
fi
if [[ ! -f "$REF_GENOME_FASTA" ]]; then
  echo "ERROR: Reference genome FASTA file '$REF_GENOME_FASTA' not found."
  exit 1
fi

# Ensure the .fai index exists for the REF_GENOME_FASTA, as mosdepth with -f will look for it
REF_GENOME_FAI="${REF_GENOME_FASTA}.fai"
if [[ ! -f "$REF_GENOME_FAI" ]]; then
  echo "INFO: FASTA index '$REF_GENOME_FAI' not found. Attempting to create it with 'samtools faidx'..."
  samtools faidx "$REF_GENOME_FASTA"
  if [[ ! -f "$REF_GENOME_FAI" ]]; then
    echo "ERROR: Failed to create FASTA index '$REF_GENOME_FAI'. Please create it manually."
    exit 1
  fi
  echo "INFO: FASTA index '$REF_GENOME_FAI' created successfully."
else
  echo "INFO: FASTA index '$REF_GENOME_FAI' already exists."
fi


mkdir -p "$COVERAGE_BASE_DIR"
if [[ ! -w "$COVERAGE_BASE_DIR" ]]; then
  echo "ERROR: Cannot write to coverage base directory '$COVERAGE_BASE_DIR'"; exit 1
fi

echo "===== SCRIPT CONFIGURATION (run_mosdepth_coverage_alt.sh) ====="
echo "Aligned BAM Dir          : $ALIGNDIR"
echo "Coverage Stats Base Dir  : $COVERAGE_BASE_DIR"
echo "Reference FASTA          : $REF_GENOME_FASTA (will look for .fai)"
echo "Threads                  : $THREADS"
echo "========================================================"

#### MOSDEPTH LOOP ####
for SAMPLE in wt1 wt2 wt3 wt4 mt1 mt2 mt3 mt4; do
  echo "--- Calculating coverage for sample $SAMPLE ---"

  INPUT_BAM="$ALIGNDIR/${SAMPLE}.aligned.sorted.bam"
  OUTPUT_PREFIX_PATH="${COVERAGE_BASE_DIR}/${SAMPLE}"

  if [[ ! -f "$INPUT_BAM" ]]; then
    echo "WARNING: Aligned BAM file '$INPUT_BAM' not found for sample $SAMPLE. Skipping."
    continue
  fi
  if [[ ! -f "${INPUT_BAM}.bai" && ! -f "${INPUT_BAM%.bam}.bai" ]]; then
    echo "WARNING: Index file for '$INPUT_BAM' not found. Please index it first. Skipping $SAMPLE."
    continue
  fi

  echo "Input BAM      : $INPUT_BAM"
  echo "Output Prefix  : $OUTPUT_PREFIX_PATH"

  # --- Attempt 1: Using mosdepth with -f FASTA_FILE ---
  # mosdepth with -f will look for a .fai index in the same location as the FASTA file.
  echo "Running mosdepth for $SAMPLE (using -f $REF_GENOME_FASTA)..."
  mosdepth -t "$THREADS" -f "$REF_GENOME_FASTA" "$OUTPUT_PREFIX_PATH" "$INPUT_BAM"

  if [ $? -eq 0 ]; then
    echo "SUCCESS: mosdepth completed for $SAMPLE using -f. Coverage reports saved with prefix $OUTPUT_PREFIX_PATH"
  else
    echo "ERROR: mosdepth failed for $SAMPLE using -f. (See error above or in nohup.err if used)"
    echo "This might still indicate an issue with the .fai file or the FASTA itself."
    continue # Skip to next sample if the primary method fails
  fi

  echo "--- Finished coverage calculation for sample $SAMPLE ---"
done

echo "========================================================"
echo "Mosdepth coverage calculation finished for all samples."
echo "Coverage reports are located in subdirectories within: $COVERAGE_BASE_DIR"
echo "Key summary file for each sample: e.g., ${COVERAGE_BASE_DIR}/wt1.mosdepth.summary.txt"
echo "========================================================"
echo "INFO: Script finished at $(date)"
