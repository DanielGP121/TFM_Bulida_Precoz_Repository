#!/usr/bin/env bash
# FINAL SCRIPT: Prepares all data tables for the GATK analysis notebook.
# It correctly extracts VAF data from the full, pre-filtered VCF.

set -euo pipefail

## CONFIGURATION ##
EXP_DIR="/data/training2/analisis_TFM_Bulida_Precoz/08_parameter_change_analysis"
GATK_WORKFLOW_DIR="${EXP_DIR}/GATK_HaplotypeCaller_workflow"

# Input VCF with the 224 strict candidates
STRICT_VCF="${GATK_WORKFLOW_DIR}/03_strict_candidate_search/strict_candidates_15pct.vcf.gz"
# Input VCF with all samples, filtered by the 15% proportional threshold
FULL_VCF="${GATK_WORKFLOW_DIR}/01_proportional_filtering/final_bcftools_propDP_15pct.filtered.vcf.gz"

# Output directory for the notebook data
NOTEBOOK_DATA_DIR="${GATK_WORKFLOW_DIR}/notebook_data"
mkdir -p "${NOTEBOOK_DATA_DIR}"

SNPEFF_GENOME_NAME="bulida_v1"

# --- Output Files ---
ANNOTATED_VCF="${NOTEBOOK_DATA_DIR}/gatk_15pct_strict.ann.vcf"
FULL_TABLE_OUTPUT="${NOTEBOOK_DATA_DIR}/gatk_15pct_strict_full_table.tsv"
VAF_TABLE_OUTPUT="${NOTEBOOK_DATA_DIR}/gatk_15pct_vaf_data.tsv"

echo "--- Preparing all data files for final GATK Notebook Analysis ---"

# --- 1. Annotate VCF with SnpEff ---
echo "Step 1: Annotating the strict candidate VCF..."
snpEff ann -v "${SNPEFF_GENOME_NAME}" "${STRICT_VCF}" > "${ANNOTATED_VCF}"

# --- 2. Extract Full Annotation Table ---
echo "Step 2: Extracting full annotation table..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' "${ANNOTATED_VCF}" > "${FULL_TABLE_OUTPUT}"

# --- 3. Extract VAF Data Table from the correct source ---
echo "Step 3: Extracting Allele Depth data for VAF analysis..."
# The -R flag tells bcftools to only process regions listed in the strict VCF file,
# but it extracts the data from the FULL VCF that contains all 8 samples.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP]\n' -R "${STRICT_VCF}" "${FULL_VCF}" > "${VAF_TABLE_OUTPUT}"
echo "VAF data table created successfully."

# --- 4. Cleanup ---
rm "${ANNOTATED_VCF}"

echo -e "\n### Data preparation complete. ###"