#!/usr/bin/env bash
# FINAL SCRIPT: Performs iterative filtering by piping a series of simple bcftools commands.
# This is the most robust, albeit less elegant, method to avoid all shell parsing errors.

set -euo pipefail

## CONFIGURATION ##
EXP_DIR="/data/training2/analisis_TFM_Bulida_Precoz/08_parameter_change_analysis"
INPUT_VCF="${EXP_DIR}/01_proportional_filtering/all_samples_8.joint.vcf.gz"
COVERAGE_SUMMARY_FILE="${EXP_DIR}/all_samples.mosdepth_summary_combined.txt"
OUTPUT_DIR="${EXP_DIR}/01_proportional_filtering"

PERCENTAGES=(3 5 7 10 15)
VCF_SAMPLE_ORDER=("mt1" "mt2" "mt3" "mt4" "wt1" "wt2" "wt3" "wt4")
THREADS=8

echo "--- Starting Final, Step-by-Step Iterative Filtering Workflow ---"

# --- 1. Extract mean coverages ---
declare -A COVERAGES
while read -r sample coverage; do
    COVERAGES[$sample]=$coverage
done < <(awk '/^=== Content of: / {split($4, a, "."); sample=a[1]} /^[t]otal/ {print sample, $4}' "${COVERAGE_SUMMARY_FILE}")

# --- 2. Main loop to process each percentage ---
for percent in "${PERCENTAGES[@]}"; do
    echo -e "\n======================================================="
    echo "--- Processing threshold: ${percent}% ---"
    
    OUTPUT_VCF_FILTERED="${OUTPUT_DIR}/final_bcftools_propDP_${percent}pct.filtered.vcf.gz"

    if [ -f "${OUTPUT_VCF_FILTERED}" ]; then
        echo "Output file for ${percent}% already exists. Skipping."
        continue
    fi

    # Calculate DP thresholds
    DP_THRESHOLDS=()
    for sample in "${VCF_SAMPLE_ORDER[@]}"; do
        coverage=${COVERAGES[$sample]}
        threshold=$(LC_NUMERIC=C printf "%.0f" $(echo "${coverage} * ${percent} / 100" | bc -l))
        DP_THRESHOLDS+=($threshold)
    done
    echo "Calculated DP thresholds for ${percent}%: ${DP_THRESHOLDS[*]}"

    # --- 3. Build and execute the command chain ---
    # Start with the input file
    CMD="bcftools view ${INPUT_VCF}"

    # Add a separate filter for EACH sample's DP
    for i in "${!VCF_SAMPLE_ORDER[@]}"; do
        CMD+=" | bcftools filter --threads ${THREADS} -i 'FORMAT/DP[${i}] >= ${DP_THRESHOLDS[${i}]}'"
    done
    
    # Add the final quality filter and redirect to output file
    CMD+=" | bcftools filter --threads ${THREADS} -i 'QD >= 2.0 && MQ >= 40.0 && FS <= 60.0' -Oz -o ${OUTPUT_VCF_FILTERED}"
    
    echo "Applying chained filters..."
    # Use eval to execute the command string with all its pipes
    eval ${CMD}

    bcftools index --threads "${THREADS}" "${OUTPUT_VCF_FILTERED}"
    
    COUNT=$(bcftools view -H "${OUTPUT_VCF_FILTERED}" | wc -l)
    echo "Filtering for ${percent}% complete. Found ${COUNT} variants."
    echo "Output file: ${OUTPUT_VCF_FILTERED}"
done

echo -e "\n### Iterative filtering workflow finished. ###"