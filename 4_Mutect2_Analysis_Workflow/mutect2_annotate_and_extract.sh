#!/usr/bin/env bash
# Annotate and extract the final list of independent Mutect2 candidates.

set -euo pipefail

## CONFIGURATION ##
BASE_DIR="/data/training2/analisis_TFM_Bulida_Precoz"
EXP_DIR="${BASE_DIR}/08_parameter_change_analysis"
# --- Input/Output Directories for the independent Mutect2 analysis ---
CANDIDATE_DIR="${EXP_DIR}/mutect2_independent_analysis/candidates"
FINAL_RESULTS_DIR="${EXP_DIR}/mutect2_independent_analysis/final_results"
mkdir -p "${FINAL_RESULTS_DIR}"

# --- SnpEff Database Name ---
SNPEFF_GENOME_NAME="bulida_v1"

# --- Input File ---
# This is the VCF containing the strict 4-out-of-4 candidates from the M2 script
INPUT_VCF="${CANDIDATE_DIR}/mutect2_strict_4of4_candidates.vcf.gz"

# --- Output Files ---
ANNOTATED_VCF="${FINAL_RESULTS_DIR}/mutect2_final.ann.vcf"
TABLE_OUTPUT="${FINAL_RESULTS_DIR}/mutect2_final_table.tsv"
GENE_LIST_OUTPUT="${FINAL_RESULTS_DIR}/mutect2_final_gene_list.txt"


echo "--- Final Annotation and Extraction for Mutect2 Candidates ---"

if [ ! -f "${INPUT_VCF}" ]; then
    echo "ERROR: Input file ${INPUT_VCF} not found. Please run script M2 first."
    exit 1
fi

# Step 1: Functional Annotation with SnpEff
echo "Step 1: Annotating final candidates..."
snpEff ann -v "${SNPEFF_GENOME_NAME}" "${INPUT_VCF}" > "${ANNOTATED_VCF}"
echo "Annotation complete."

# Step 2: Extraction to a TSV table
# We extract all candidates, as they have already passed the strict biological filter.
echo "Step 2: Extracting all candidates into a table..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' "${ANNOTATED_VCF}" > "${TABLE_OUTPUT}"
echo "Final Mutect2 table created: ${TABLE_OUTPUT}"

# Step 3: Create final gene list
echo "Step 3: Creating final gene list..."
cut -f 5 "${TABLE_OUTPUT}" | cut -d '|' -f 5 | sort -u | grep -v '^\s*$' > "${GENE_LIST_OUTPUT}"
echo "Final Mutect2 gene list created: ${GENE_LIST_OUTPUT}"

echo -e "\n### Independent Mutect2 analysis workflow finished! ###"
echo "Final results are located in: ${FINAL_RESULTS_DIR}"