#!/usr/bin/env bash
# SCRIPT M2: Applies a proportional DP filter and keeps only PASSing variants for the mutant sample.

set -euo pipefail

## CONFIGURATION ##
# --- !! SET THE DESIRED PERCENTAGE THRESHOLD HERE !! ---
PERCENT_TO_USE=15 # Using 15% as a robust, conservative choice

BASE_DIR="/data/training2/analisis_TFM_Bulida_Precoz"
EXP_DIR="${BASE_DIR}/08_parameter_change_analysis"
# --- Corrected input directory name ---
INPUT_DIR="${EXP_DIR}/mutect2_final_workflow/01_mutect2_filtered"
COVERAGE_SUMMARY_FILE="${EXP_DIR}/all_samples.mosdepth_summary_combined.txt"
OUTPUT_DIR="${EXP_DIR}/mutect2_final_workflow/02_dp_filtered"
mkdir -p "${OUTPUT_DIR}"

MUTANT_SAMPLES=("mt1" "mt2" "mt3" "mt4")
THREADS=8

# --- 1. Calculate DP thresholds ---
declare -A DP_THRESHOLDS
while read -r sample coverage; do
    threshold=$(LC_NUMERIC=C printf "%.0f" $(echo "${coverage} * ${PERCENT_TO_USE} / 100" | bc -l))
    DP_THRESHOLDS[$sample]=$threshold
done < <(awk '/^=== Content of: / {split($4, a, "."); sample=a[1]} /^[t]otal/ {print sample, $4}' "${COVERAGE_SUMMARY_FILE}")

# --- 2. Main loop to filter each sample ---
echo "--- Applying Proportional DP Filter (${PERCENT_TO_USE}%) ---"
for sample in "${MUTANT_SAMPLES[@]}"; do
    INPUT_VCF="${INPUT_DIR}/${sample}.filtered.vcf.gz"
    OUTPUT_VCF="${OUTPUT_DIR}/${sample}.propDP_${PERCENT_TO_USE}pct.vcf.gz"
    DP_THRESHOLD=${DP_THRESHOLDS[$sample]}

    if [ -f "${OUTPUT_VCF}" ]; then
        echo "Final filtered file for ${sample} already exists. Skipping."
        continue
    fi

    echo "Processing ${sample} with DP threshold >= ${DP_THRESHOLD}..."
    
    # This command chain does everything:
    # 1. Selects only the mutant sample column.
    # 2. Keeps only variants marked as PASS from the previous GATK step.
    # 3. Keeps only variants where that sample's DP meets the proportional threshold.
    bcftools view --threads "${THREADS}" --samples "${sample}" "${INPUT_VCF}" \
    | bcftools filter --threads "${THREADS}" -i 'FILTER="PASS"' \
    | bcftools filter --threads "${THREADS}" -i "FORMAT/DP >= ${DP_THRESHOLD}" -Oz -o "${OUTPUT_VCF}"

    bcftools index --threads "${THREADS}" "${OUTPUT_VCF}"
done
echo "Proportional DP filtering complete. Final VCFs are in ${OUTPUT_DIR}"