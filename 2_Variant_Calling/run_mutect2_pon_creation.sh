#!/usr/bin/env bash

# run_mutect2_pon_creation.sh

# Create a Panel of Normals (PoN) by merging filtered single-sample VCFs
# using bcftools, then calling CreateSomaticPanelOfNormals on the merged VCF.

set -euo pipefail
trap 'echo "ERROR on line $LINENO: exit code $?" >&2; exit 1' ERR

# --- CONFIGURATION ---
# Directory containing deduplicated BAMs
DEDUP_DIR="/data/training2/analisis_TFM_Bulida_Precoz/03_deduplication_picard"

# Main output directory for Mutect2 and PoN
PON_DIR="/data/training2/analisis_TFM_Bulida_Precoz/04_mutect2_pon"

# Reference FASTA
REF_GENOME="/data/training2/info/assemblies/BUL_cur_guided.v1.0.fasta"

# Full path to GATK jar
GATK_JAR="/data/training2/softwares/conda_envs/env_tfm_bulida/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"

# List of normal sample IDs (no extensions)
NORMAL_SAMPLES="wt1 wt2 wt3 wt4"

# JVM memory (MB) and threads for pair-HMM
MEM_MB=20480
THREADS=8

# --- PREPARE DIRECTORIES ---
mkdir -p \
  "${PON_DIR}/single_sample_vcfs" \
  "${PON_DIR}/normalized_vcfs" \
  "${PON_DIR}/filtered_vcfs"

# --- STEP 1: Run Mutect2 in tumor-only mode on each normal sample ---
echo "=== STEP 1: Mutect2 (tumor-only) on normals ==="
for SAMPLE in ${NORMAL_SAMPLES}; do
  BAM="${DEDUP_DIR}/${SAMPLE}.dedup.bam"
  RAW_VCF="${PON_DIR}/single_sample_vcfs/${SAMPLE}.single_sample.vcf.gz"
  if [[ ! -f "${BAM}" ]]; then
    echo "WARNING: BAM not found (${BAM}), skipping ${SAMPLE}"
    continue
  fi
  if [[ -f "${RAW_VCF}" ]]; then
    echo "INFO: VCF already exists (${RAW_VCF}), skipping"
    continue
  fi

  gatk --java-options "-Xmx${MEM_MB}m" Mutect2 \
    -R "${REF_GENOME}" \
    -I "${BAM}" \
    -O "${RAW_VCF}" \
    --native-pair-hmm-threads ${THREADS}
done

# --- STEP 2: Normalize VCFs and split MNPs ---
echo -e "\n=== STEP 2: Normalize & split MNPs ==="
for SAMPLE in ${NORMAL_SAMPLES}; do
  RAW_VCF="${PON_DIR}/single_sample_vcfs/${SAMPLE}.single_sample.vcf.gz"
  NORM_VCF="${PON_DIR}/normalized_vcfs/${SAMPLE}.normalized.vcf.gz"
  if [[ -f "${NORM_VCF}" ]]; then
    echo "INFO: Normalized VCF exists (${NORM_VCF}), skipping"
    continue
  fi

  gatk --java-options "-Xmx${MEM_MB}m" LeftAlignAndTrimVariants \
    -R "${REF_GENOME}" \
    -V "${RAW_VCF}" \
    -O "${NORM_VCF}" \
    --split-multi-allelics \
    --split-mnp

  gatk IndexFeatureFile -I "${NORM_VCF}"
done

# --- STEP 3: Filter out MNPs (keep only SNPs and INDELs) ---
echo -e "\n=== STEP 3: Filter out MNPs ==="
for SAMPLE in ${NORMAL_SAMPLES}; do
  NORM_VCF="${PON_DIR}/normalized_vcfs/${SAMPLE}.normalized.vcf.gz"
  FILT_VCF="${PON_DIR}/filtered_vcfs/${SAMPLE}.filtered.vcf.gz"
  if [[ -f "${FILT_VCF}" ]]; then
    echo "INFO: Filtered VCF exists (${FILT_VCF}), skipping"
    continue
  fi

  gatk --java-options "-Xmx${MEM_MB}m" SelectVariants \
    -R "${REF_GENOME}" \
    -V "${NORM_VCF}" \
    -O "${FILT_VCF}" \
    --select-type-to-include SNP \
    --select-type-to-include INDEL

  gatk IndexFeatureFile -I "${FILT_VCF}"
done

# --- STEP 4: Merge all filtered VCFs into a single multi-sample VCF using bcftools ---
echo -e "\n=== STEP 4: bcftools merge ==="
MERGED_VCF="${PON_DIR}/all_normals_for_pon.vcf.gz"

if [[ ! -f "${MERGED_VCF}" ]]; then
  bcftools merge \
    $(for SAMPLE in ${NORMAL_SAMPLES}; do
        echo -n "${PON_DIR}/filtered_vcfs/${SAMPLE}.filtered.vcf.gz "
      done) \
    -Oz -o "${MERGED_VCF}" \
    --threads ${THREADS}

  tabix -p vcf "${MERGED_VCF}"
fi

# --- STEP 5: Create the Panel of Normals from the merged VCF ---
echo -e "\n=== STEP 5: CreateSomaticPanelOfNormals ==="
gatk --java-options "-Xmx${MEM_MB}m" CreateSomaticPanelOfNormals \
  -R "${REF_GENOME}" \
  -V "${MERGED_VCF}" \
  -O "${PON_DIR}/bulida_pon.vcf.gz"

echo -e "\nPanel of Normals successfully created at:\n${PON_DIR}/bulida_pon.vcf.gz"







