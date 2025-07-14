#!/usr/bin/env bash
# SCRIPT: Iterates through proportionally filtered VCFs to find candidates
# that match the strict biological criteria (4/4 in mutants, 0/4 in wild-types).

set -euo pipefail

## CONFIGURATION ##
EXP_DIR="/data/training2/analisis_TFM_Bulida_Precoz/08_parameter_change_analysis"
INPUT_DIR="${EXP_DIR}/01_proportional_filtering"
# A new directory to store the final results of this analysis
OUTPUT_DIR="${EXP_DIR}/03_strict_candidate_search"
mkdir -p "${OUTPUT_DIR}"

PERCENTAGES=(3 5 7 10 15)
VCF_SAMPLE_ORDER=("mt1" "mt2" "mt3" "mt4" "wt1" "wt2" "wt3" "wt4")
THREADS=8

# --- Define the strict biological filter expression ---
# Include variants where ALL 4 mutants have the variant AND ALL 4 wild-types do NOT.
STRICT_FILTER_EXP='N_PASS(GT[0-3]="alt")==4 && N_PASS(GT[4-7]="alt")==0'

## SCRIPT LOGIC ##
echo "--- Starting Workflow: Searching for Strict 4/4 Candidates ---"
echo "Filter to be applied: ${STRICT_FILTER_EXP}"

# Create a summary file header
SUMMARY_FILE="${OUTPUT_DIR}/strict_candidate_summary.txt"
echo -e "Threshold_Percent\tFinal_Candidate_Count" > "${SUMMARY_FILE}"


for percent in "${PERCENTAGES[@]}"; do
    echo -e "\n======================================================="
    echo "--- Analyzing file from threshold: ${percent}% ---"
    
    INPUT_VCF="${INPUT_DIR}/final_bcftools_propDP_${percent}pct.filtered.vcf.gz"
    OUTPUT_VCF="${OUTPUT_DIR}/strict_candidates_${percent}pct.vcf.gz"

    if [ ! -f "${INPUT_VCF}" ]; then
        echo "Input file ${INPUT_VCF} not found. Skipping."
        continue
    fi

    # Apply the strict biological filter
    bcftools view \
      --threads "${THREADS}" \
      --samples "${VCF_SAMPLE_ORDER}" \
      --include "${STRICT_FILTER_EXP}" \
      -Oz -o "${OUTPUT_VCF}" \
      "${INPUT_VCF}"
    
    # Count the results
    COUNT=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
    
    echo "Result for ${percent}% threshold: Found ${COUNT} strict candidates."
    
    # Append the result to the summary file
    echo -e "${percent}\t${COUNT}" >> "${SUMMARY_FILE}"
done

echo -e "\n### Workflow Complete! ###"
echo "Summary of results:"
cat "${SUMMARY_FILE}"