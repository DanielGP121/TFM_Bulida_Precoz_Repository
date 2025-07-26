#!/usr/bin/env bash
# SCRIPT: Finds variants present in all wild-type samples (4/4) and
# absent from all mutant samples (0/4), and then performs functional annotation.
# This version saves all outputs to the current working directory.

set -euo pipefail

## ============================================================================================
## CONFIGURATION
## ============================================================================================

# --- Input File ---
# The VCF that has passed the rigorous 15% proportional depth filter and contains all 8 samples.
INPUT_VCF="/data/training2/analisis_TFM_Bulida_Precoz/08_parameter_change_analysis/GATK_HaplotypeCaller_workflow/01_proportional_filtering/final_bcftools_propDP_15pct.filtered.vcf.gz"

# --- Output Directory ---
# Use the current directory for all outputs.
OUTPUT_DIR=$(pwd)

# --- Parameters ---
# The order must match the VCF header: mt1,mt2,mt3,mt4,wt1,wt2,wt3,wt4
VCF_SAMPLE_ORDER="mt1,mt2,mt3,mt4,wt1,wt2,wt3,wt4"
# This filter selects variants present in all 4 wild-types (indices 4-7)
# and homozygous reference in all 4 mutants (indices 0-3).
STRICT_FILTER_EXP_WT='N_PASS(GT[4-7]="alt")==4 && N_PASS(GT[0-3]="ref")==4'
THREADS=8
SNPEFF_GENOME_NAME="bulida_v1"

# --- Output Files (named for the current directory) ---
OUTPUT_VCF_WT="${OUTPUT_DIR}/inverse_candidates_0_vs_4.vcf.gz"
TABLE_OUTPUT_WT="${OUTPUT_DIR}/inverse_candidates_0_vs_4_annotated_table.tsv"
GENE_LIST_OUTPUT_WT="${OUTPUT_DIR}/inverse_candidates_0_vs_4_gene_list.txt"


## ============================================================================================
## SCRIPT LOGIC
## ============================================================================================

echo "--- Starting Workflow: Finding Wild-Type Specific (Inverse) Variants ---"

if [ ! -f "${INPUT_VCF}" ]; then
    echo "ERROR: Input VCF not found at ${INPUT_VCF}"
    exit 1
fi

# --- Step 1: Apply the strict biological filter for wild-type variants ---
echo "Step 1: Applying biological filter: ${STRICT_FILTER_EXP_WT}"
bcftools view \
    --threads "${THREADS}" \
    --samples "${VCF_SAMPLE_ORDER}" \
    --include "${STRICT_FILTER_EXP_WT}" \
    -Oz -o "${OUTPUT_VCF_WT}" \
    "${INPUT_VCF}"

bcftools index --threads "${THREADS}" "${OUTPUT_VCF_WT}"
COUNT=$(bcftools view -H "${OUTPUT_VCF_WT}" | wc -l)
echo "Found ${COUNT} variants specific to the wild-type group."

# --- Step 2: Functional Annotation with SnpEff ---
if [ "${COUNT}" -gt 0 ]; then
    echo -e "\nStep 2: Annotating variants with SnpEff..."
    # Create a temporary annotated VCF file
    ANNOTATED_VCF_WT_TEMP="${OUTPUT_DIR}/temp_inverse_annotated.vcf"
    snpEff ann -v "${SNPEFF_GENOME_NAME}" "${OUTPUT_VCF_WT}" > "${ANNOTATED_VCF_WT_TEMP}"

    # --- Step 3: Extract results to table and gene list ---
    echo "Step 3: Extracting results to table and gene list..."
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' "${ANNOTATED_VCF_WT_TEMP}" > "${TABLE_OUTPUT_WT}"
    cut -f 5 "${TABLE_OUTPUT_WT}" | cut -d '|' -f 5 | sort -u | grep -v '^\s*$' > "${GENE_LIST_OUTPUT_WT}"
    
    # Cleanup temporary file
    rm "${ANNOTATED_VCF_WT_TEMP}"
else
    echo "No wild-type specific variants found. Skipping annotation."
fi

echo -e "\n### Workflow Complete! ###"
echo "Final results are in this directory: ${OUTPUT_DIR}"