#!/usr/bin/env bash
set -euo pipefail
trap 'echo "ERROR: Script failed on line $LINENO with exit code $?"; exit 1' ERR

# run_mutect2_pooled_wt.sh
# Merges all wild-type (wt) BAMs to create a pooled control sample,
# then runs GATK Mutect2 for each mutant sample against this pooled wild-type.
# CORRECTED to pass all normal sample names to Mutect2.

# ── CONFIGURATION ────────────────────────────────────────────────────────────────
DEDUPDIR="/data/training2/analisis_TFM_Bulida_Precoz/03_deduplication_picard" 
MUTECT2_DIR="/data/training2/analisis_TFM_Bulida_Precoz/04_mutect2_calling" 
REF_GENOME="/data/training2/info/assemblies/BUL_cur_guided.v1.0.fasta" 
PON_VCF="/data/training2/analisis_TFM_Bulida_Precoz/04_mutect2_pon/bulida_pon.vcf.gz" 

# GATK and Java options
MEM_MB=20480   # Approx 20GB [cite: User's Request]
THREADS=8      # Used for samtools merge [cite: User's Request]

# Wild-type and mutant sample lists
WILD_TYPE_SAMPLES="wt1 wt2 wt3 wt4" 
MUTANT_SAMPLES="mt1 mt2 mt3 mt4" 

# Define the name for the pooled wild-type BAM file
POOLED_WT_BAM="${DEDUPDIR}/pooled_wild_type.dedup.bam" 

# ── PREPARATION ───────────────────────────────────────────────────────────────────
echo "INFO: GATK Mutect2 Pooled Wild-Type script started at $(date)" 

# ... (Prerequisite checks for gatk, samtools, and directories would be here) ...

# ── STAGE 1: MERGE WILD-TYPE BAMS ────────────────────────────────────────
echo "--- Stage 1: Merging wild-type BAMs to create a pooled control sample ---" 

WT_BAM_LIST=()
for SAMPLE in $WILD_TYPE_SAMPLES; do 
    WT_BAM_PATH="${DEDUPDIR}/${SAMPLE}.dedup.bam" 
    if [[ -f "$WT_BAM_PATH" ]]; then 
        WT_BAM_LIST+=("$WT_BAM_PATH") 
    else
        echo "WARNING: Wild-type BAM not found for ${SAMPLE}, excluding." ]
    fi
done

if [ ${#WT_BAM_LIST[@]} -eq 0 ]; then 
    echo "ERROR: No wild-type BAM files found to create a pool." #
    exit 1 
fi

echo "  Merging into ${POOLED_WT_BAM}:" 
printf "    %s\n" "${WT_BAM_LIST[@]}" 

# Use samtools merge: output BAM first, then inputs. -r attaches RG from first file. -h includes headers from all.
(samtools merge -@ $THREADS -r -h "${WT_BAM_LIST[0]}" -f "$POOLED_WT_BAM" "${WT_BAM_LIST[@]}" && \
  echo "  SUCCESS: Merged wild-type BAMs." && \
  samtools index "$POOLED_WT_BAM" && \
  echo "  SUCCESS: Pooled wild-type BAM indexed." \
) || { echo "  ERROR: Failed to merge or index wild-type BAMs."; exit 1; } 

# ── STAGE 2: MUTECT2 VARIANT CALLING & FILTERING ──────────────────────────────────────────
echo "--- Stage 2: Running Mutect2 for each mutant against the pooled wild-type ---" 
for MUTANT_SAMPLE in $MUTANT_SAMPLES; do 
  echo
  echo "--- Processing mutant sample: $MUTANT_SAMPLE vs Pooled Wild-Type ---" 

  MUTANT_BAM="${DEDUPDIR}/${MUTANT_SAMPLE}.dedup.bam" 
  
  RAW_VCF_OUT="${MUTECT2_DIR}/raw_vcfs/${MUTANT_SAMPLE}_vs_pooled_wt.raw.vcf.gz" 
  FILTERED_VCF_OUT="${MUTECT2_DIR}/filtered_vcfs/${MUTANT_SAMPLE}_vs_pooled_wt.filtered.vcf.gz" 

  [[ -f "$MUTANT_BAM" ]] || { echo "WARNING: $MUTANT_BAM not found, skipping."; continue; } 

  # Build an array of normal sample arguments for Mutect2
  NORMAL_ARGS=()
  for WT_SAMPLE in $WILD_TYPE_SAMPLES; do
      NORMAL_ARGS+=("-normal" "$WT_SAMPLE")
  done

  # Step 1: Run Mutect2
  echo "  Step 1/3: Running Mutect2 for ${MUTANT_SAMPLE}..." 
  (gatk --java-options "-Xmx${MEM_MB}m" Mutect2 \
    -R "$REF_GENOME" \
    -I "$MUTANT_BAM" \
    -I "$POOLED_WT_BAM" \
    "${NORMAL_ARGS[@]}" \
    -pon "$PON_VCF" \
    --f1r2-tar-gz "${MUTECT2_DIR}/stats/${MUTANT_SAMPLE}.f1r2.tar.gz" \
    -O "$RAW_VCF_OUT" && \
  echo "  SUCCESS: Mutect2 completed for ${MUTANT_SAMPLE}" \
  ) || { echo "  ERROR: Mutect2 failed for ${MUTANT_SAMPLE}."; continue; } 

  # Step 2: GetPileupSummaries
  echo "  Step 2/3: Running GetPileupSummaries for ${MUTANT_SAMPLE}..." 
  (gatk GetPileupSummaries \
    -I "$MUTANT_BAM" \
    -V "$RAW_VCF_OUT" \
    -L "$RAW_VCF_OUT" \
    -O "${MUTECT2_DIR}/stats/${MUTANT_SAMPLE}.pileups.table" && \
  echo "  SUCCESS: GetPileupSummaries completed for ${MUTANT_SAMPLE}." \
  ) || echo "  WARNING: GetPileupSummaries failed; contamination estimate may be off." 

  # Step 3: FilterMutectCalls
  echo "  Step 3/3: Running FilterMutectCalls for ${MUTANT_SAMPLE}..." 
  (gatk FilterMutectCalls \
    -V "$RAW_VCF_OUT" \
    --contamination-table "${MUTECT2_DIR}/stats/${MUTANT_SAMPLE}.pileups.table" \
    -O "$FILTERED_VCF_OUT" && \
  echo "  SUCCESS: FilterMutectCalls completed for ${MUTANT_SAMPLE}" \
  ) || { echo "  ERROR: FilterMutectCalls failed for ${MUTANT_SAMPLE}."; continue; } 

done

echo "============================================================" 
echo "Mutect2 calling and filtering finished for all mutant samples." 
echo "Filtered VCFs in: ${MUTECT2_DIR}/filtered_vcfs/" 
echo "Script finished at $(date)" 