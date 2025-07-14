#!/usr/bin/env bash
# Applies the default GATK FilterMutectCalls filters to the raw VCFs.

set -euo pipefail

## CONFIGURATION ##
BASE_DIR="/data/training2/analisis_TFM_Bulida_Precoz"
EXP_DIR="${BASE_DIR}/08_parameter_change_analysis"
RAW_VCF_DIR="${BASE_DIR}/04_mutect2_calling/raw_vcfs"
# A new directory for the outputs of this workflow
OUTPUT_DIR="${EXP_DIR}/mutect2_final_workflow"

# --- CORRECTED DIRECTORY NAME AS PER YOUR SUGGESTION ---
FILTERED_DIR="${OUTPUT_DIR}/01_mutect2_filtered"
mkdir -p "${FILTERED_DIR}"

REF_GENOME="/data/training2/info/assemblies/BUL_cur_guided.v1.0.fasta"
MUTANT_SAMPLES=("mt1" "mt2" "mt3" "mt4")

echo "--- Running GATK FilterMutectCalls on Raw VCFs ---"
for sample in "${MUTANT_SAMPLES[@]}"; do
    RAW_VCF="${RAW_VCF_DIR}/${sample}_vs_pooled_wt.raw.vcf.gz"
    FILTERED_VCF_OUT="${FILTERED_DIR}/${sample}.filtered.vcf.gz"
    
    if [ -f "${FILTERED_VCF_OUT}" ]; then
        echo "Filtered file for ${sample} already exists. Skipping."
        continue
    fi
    
    echo "Processing ${sample}..."
    gatk FilterMutectCalls -R "${REF_GENOME}" -V "${RAW_VCF}" -O "${FILTERED_VCF_OUT}"
done
echo "GATK filtering complete. Tagged VCFs are in: ${FILTERED_DIR}"