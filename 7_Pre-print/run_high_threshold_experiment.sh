#!/usr/bin/env bash
# SCRIPT: Performs an iterative filtering experiment with an expanded range of DP thresholds (3-50%) on the Mutect2 results.
set -euo pipefail

## ==================
## CONFIGURATION
## ==================

# --- Thresholds to test ---
PERCENTAGES=(3 5 7 12 15 20 25 30 35 50)

# --- Paths to INPUT files (based on your file hierarchy) ---
COVERAGE_SUMMARY_FILE="/data/training2/analisis_TFM_Bulida_Precoz/08_parameter_change_analysis/all_samples.mosdepth_summary_combined.txt"
MUTECT2_FILTERED_DIR="/data/training2/analisis_TFM_Bulida_Precoz/08_parameter_change_analysis/mutect2_final_workflow/01_mutect2_filtered"

# --- OUTPUT directories (inside the current '7_Pre-print' directory) ---
OUTPUT_BASE_DIR=$(pwd) # Use the current directory as the base
DP_FILTERED_DIR="${OUTPUT_BASE_DIR}/01_dp_filtered_expanded_thresholds"
FINAL_CANDIDATES_DIR="${OUTPUT_BASE_DIR}/02_final_candidates_expanded_thresholds"

mkdir -p "${DP_FILTERED_DIR}"
mkdir -p "${FINAL_CANDIDATES_DIR}"

# --- General Parameters ---
MUTANT_SAMPLES=("mt1" "mt2" "mt3" "mt4")
THREADS=8
STRICT_FILTER_EXP='N_PASS(GT="alt")==4'

echo "--- Starting expanded-threshold filtering experiment ---"

## ==================
## STEP A: Calculate Depth (DP) thresholds for each sample
## ==================
declare -A DP_THRESHOLDS
while read -r sample coverage; do
    DP_THRESHOLDS[$sample]=$coverage
done < <(awk '/^=== Content of: / {split($4, a, "."); sample=a[1]} /^[t]otal/ {print sample, $4}' "${COVERAGE_SUMMARY_FILE}")
echo "Mean coverages extracted for threshold calculation."


## ==================
## STEP B: Main loop to iterate through each percentage
## ==================
SUMMARY_FILE="${FINAL_CANDIDATES_DIR}/expanded_threshold_candidate_summary.txt"
echo -e "Threshold_Percent\tFinal_Candidate_Count" > "${SUMMARY_FILE}"

for percent in "${PERCENTAGES[@]}"; do
    echo -e "\n======================================================="
    echo "--- Processing threshold: ${percent}% ---"

    # Directory for VCFs from this specific threshold
    CURRENT_DP_DIR="${DP_FILTERED_DIR}/${percent}pct"
    mkdir -p "${CURRENT_DP_DIR}"

    # --- B1: Apply the proportional DP filter to each mutant sample ---
    VCF_LIST_FOR_MERGE=()
    for sample in "${MUTANT_SAMPLES[@]}"; do
        INPUT_VCF="${MUTECT2_FILTERED_DIR}/${sample}.filtered.vcf.gz"
        OUTPUT_VCF="${CURRENT_DP_DIR}/${sample}.propDP_${percent}pct.vcf.gz"
        
        # Calculate the specific DP threshold for this sample and percentage
        coverage=${DP_THRESHOLDS[$sample]}
        threshold=$(LC_NUMERIC=C printf "%.0f" $(echo "${coverage} * ${percent} / 100" | bc -l))
        
        echo "Processing ${sample} with DP threshold >= ${threshold}..."
        
        # Command chain for filtering
        bcftools view --threads "${THREADS}" --samples "${sample}" "${INPUT_VCF}" | \
        bcftools filter --threads "${THREADS}" -i 'FILTER="PASS"' | \
        bcftools filter --threads "${THREADS}" -i "FORMAT/DP >= ${threshold}" -Oz -o "${OUTPUT_VCF}"
        
        bcftools index --threads "${THREADS}" "${OUTPUT_VCF}"
        # FIX: Corrected variable name from VCF_LIST_FOR_MERge to VCF_LIST_FOR_MERGE
        VCF_LIST_FOR_MERGE+=("${OUTPUT_VCF}")
    done

    # --- B2: Merge filtered VCFs and find strict candidates ---
    STRICT_CANDIDATES_VCF="${FINAL_CANDIDATES_DIR}/strict_candidates_${percent}pct.vcf.gz"
    echo "Merging VCFs and applying biological filter (${STRICT_FILTER_EXP})..."

    # Merge and filter in a single step
    bcftools merge --threads "${THREADS}" -Oz -o - "${VCF_LIST_FOR_MERGE[@]}" | \
    bcftools view --threads "${THREADS}" --include "${STRICT_FILTER_EXP}" -Oz -o "${STRICT_CANDIDATES_VCF}"

    # --- B3: Count and report the results ---
    COUNT=$(bcftools view -H "${STRICT_CANDIDATES_VCF}" | wc -l)
    echo "Result for ${percent}% threshold: Found ${COUNT} strict candidates."
    echo -e "${percent}\t${COUNT}" >> "${SUMMARY_FILE}"
done

echo -e "\n### Filtering experiment completed successfully. ###"
echo "Summary of results:"
cat "${SUMMARY_FILE}" | column -t