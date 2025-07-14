#!/usr/bin/env bash
set -euo pipefail
trap 'echo "ERROR: Script failed on line $LINENO with exit code $?"; exit 1' ERR

# run_alignment_pipeline.sh
# Aligns trimmed paired-end FASTQ files using BWA-MEM,
# then converts SAM to sorted BAM, fixmates, and indexes.

# ── CONFIGURATION ────────────────────────────────────────────────────────────────
TRIMDIR="/data/training2/analisis_TFM_Bulida_Precoz/01_trimming"
ALIGNDIR="/data/training2/analisis_TFM_Bulida_Precoz/02_alignment"
REF_GENOME="/data/training2/info/assemblies/BUL_cur_guided.v1.0.fasta"
THREADS=8
PLATFORM="ILLUMINA"
LIB_PREFIX="lib"
LOGDIR="${ALIGNDIR}/logs" # Directorio para logs de bwa y del script

# ── PREPARE ────────────────────────────────────────────────────────────────────
echo "INFO: Script started at $(date)"

if ! command -v bwa >/dev/null; then
  echo "ERROR: bwa command not found. Activate Conda environment or install."; exit 1
fi
if ! command -v samtools >/dev/null; then
  echo "ERROR: samtools command not found. Activate Conda environment or install."; exit 1
fi

if [[ ! -d "$TRIMDIR" ]]; then
  echo "ERROR: Trimmed reads directory '$TRIMDIR' does not exist"; exit 1
fi
if [[ ! -f "$REF_GENOME" ]]; then
  echo "ERROR: Reference genome file '$REF_GENOME' not found"; exit 1
fi

mkdir -p "$ALIGNDIR" "$LOGDIR"
if [[ ! -w "$ALIGNDIR" || ! -w "$LOGDIR" ]]; then
  echo "ERROR: Cannot write to output directory '$ALIGNDIR' or log directory '$LOGDIR'"; exit 1
fi

# Index the reference genome with BWA if index files are not present
if [[ ! -f "${REF_GENOME}.bwt" ]]; then
  echo "INFO: BWA index for reference genome '${REF_GENOME}' not found. Indexing now..."
  bwa index "$REF_GENOME"
  echo "INFO: BWA index created successfully."
else
  echo "INFO: BWA index for reference genome '${REF_GENOME}' already exists."
fi

echo "===== SCRIPT CONFIGURATION (run_alignment_pipeline.sh) ====="
echo "Trimmed Reads Dir   : $TRIMDIR"
echo "Alignment Output Dir: $ALIGNDIR"
echo "Logs Dir            : $LOGDIR"
echo "Reference Genome    : $REF_GENOME"
echo "Threads             : $THREADS"
echo "======================================================"

#### ALIGNMENT AND BAM PROCESSING LOOP ####
for SAMPLE in wt1 wt2 wt3 wt4 mt1 mt2 mt3 mt4; do
  echo "--- Processing alignment for sample $SAMPLE ---"

  # Input trimmed FASTQ files
  R1_TRIMMED="$TRIMDIR/${SAMPLE}-trimmed-pair1.fastq.gz"
  R2_TRIMMED="$TRIMDIR/${SAMPLE}-trimmed-pair2.fastq.gz"

  # Intermediate and final file names for this sample
  SAM_OUT="$ALIGNDIR/${SAMPLE}.sam"
  BAM_INITIAL="$ALIGNDIR/${SAMPLE}.initial.bam" # Raw BAM from SAM
  BAM_SORTED_NAME="$ALIGNDIR/${SAMPLE}.sortedN.bam" # Sorted by name
  BAM_FIXMATE="$ALIGNDIR/${SAMPLE}.fixmate.bam" # After fixmate
  FINAL_BAM_FOR_DEDUP="$ALIGNDIR/${SAMPLE}.aligned.sorted.bam" # BAM listo para MarkDuplicates

  # Check if input trimmed files exist
  if [[ ! -f "$R1_TRIMMED" || ! -f "$R2_TRIMMED" ]]; then
    echo "WARNING: Trimmed FASTQ files for sample $SAMPLE ('$R1_TRIMMED', '$R2_TRIMMED') not found. Skipping."
    continue
  fi

  echo "Input R1: $R1_TRIMMED"
  echo "Input R2: $R2_TRIMMED"

  # Construct Read Group string (essential for GATK)
  RG_STRING="@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tLB:${LIB_PREFIX}_${SAMPLE}\\tPL:${PLATFORM}"
  echo "Read Group for BWA: ${RG_STRING}"

  # 1. BWA-MEM Alignment
  echo "Step 1/6: Running BWA-MEM for $SAMPLE..."
  bwa mem -t "$THREADS" -M -R "$RG_STRING" "$REF_GENOME" "$R1_TRIMMED" "$R2_TRIMMED" > "$SAM_OUT" 2> "${LOGDIR}/${SAMPLE}_bwa.err"
  if [ $? -eq 0 ]; then
    echo "SUCCESS: BWA-MEM completed."
  else
    echo "ERROR: BWA-MEM failed for $SAMPLE. Check '${LOGDIR}/${SAMPLE}_bwa.err'."
    rm -f "$SAM_OUT" # Clean up potentially incomplete SAM
    continue
  fi

  # 2. Convert SAM to BAM
  echo "Step 2/6: Converting SAM to BAM..."
  samtools view -@ "$THREADS" -S -b "$SAM_OUT" > "$BAM_INITIAL"
  if [ $? -eq 0 ]; then
    echo "SUCCESS: SAM to BAM conversion completed."
  else
    echo "ERROR: samtools view (SAM to BAM) failed."
    rm -f "$SAM_OUT" "$BAM_INITIAL"
    continue
  fi

  # 3. Sort BAM by read name (required for samtools fixmate)
  echo "Step 3/6: Sorting BAM by name..."
  samtools sort -@ "$THREADS" -n "$BAM_INITIAL" -o "$BAM_SORTED_NAME"
  if [ $? -eq 0 ]; then
    echo "SUCCESS: BAM sorted by name."
  else
    echo "ERROR: samtools sort -n (by name) failed."
    rm -f "$SAM_OUT" "$BAM_INITIAL" "$BAM_SORTED_NAME"
    continue
  fi

  # 4. Fix mate pair information and add ms tag
  echo "Step 4/6: Running samtools fixmate..."
  samtools fixmate -@ "$THREADS" -m -O bam "$BAM_SORTED_NAME" "$BAM_FIXMATE"
  if [ $? -eq 0 ]; then
    echo "SUCCESS: samtools fixmate completed."
  else
    echo "ERROR: samtools fixmate failed."
    rm -f "$SAM_OUT" "$BAM_INITIAL" "$BAM_SORTED_NAME" "$BAM_FIXMATE"
    continue
  fi

  # 5. Sort BAM by coordinate (standard for downstream analysis)
  echo "Step 5/6: Sorting BAM by coordinate..."
  samtools sort -@ "$THREADS" "$BAM_FIXMATE" -o "$FINAL_BAM_FOR_DEDUP"
  if [ $? -eq 0 ]; then
    echo "SUCCESS: BAM sorted by coordinate."
  else
    echo "ERROR: samtools sort (by coordinate) failed."
    rm -f "$SAM_OUT" "$BAM_INITIAL" "$BAM_SORTED_NAME" "$BAM_FIXMATE" "$FINAL_BAM_FOR_DEDUP"
    continue
  fi

  # 6. Index the final sorted BAM (Picard MarkDuplicates also needs the index of the input)
  echo "Step 6/6: Indexing final BAM for $SAMPLE..."
  samtools index "$FINAL_BAM_FOR_DEDUP"
  if [ $? -eq 0 ]; then
    echo "SUCCESS: Final BAM indexed: ${FINAL_BAM_FOR_DEDUP}.bai"
  else
    echo "ERROR: samtools index failed."
    # Consider removing FINAL_BAM_FOR_DEDUP if indexing fails and it's critical
    rm -f "$SAM_OUT" "$BAM_INITIAL" "$BAM_SORTED_NAME" "$BAM_FIXMATE" "$FINAL_BAM_FOR_DEDUP"
    continue
  fi

  # Cleanup intermediate files for this sample
  echo "Cleaning up intermediate files for $SAMPLE..."
  rm -f "$SAM_OUT" "$BAM_INITIAL" "$BAM_SORTED_NAME" "$BAM_FIXMATE"
  echo "SUCCESS: Intermediate files cleaned up."

  echo "--- Finished alignment and BAM pre-processing for sample $SAMPLE ---"
done

echo "================================================================="
echo "Alignment and BAM pre-processing finished for all samples."
echo "Final BAM files (ready for deduplication) are in: $ALIGNDIR"
echo "BWA error logs are in: $LOGDIR"
echo "================================================================="
echo "INFO: Script finished at $(date)"
