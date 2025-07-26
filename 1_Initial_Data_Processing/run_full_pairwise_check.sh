#!/usr/bin/env bash
# SCRIPT: Performs a pairwise comparison of ALL samples (mutant and wild-type).
# Corrected Version: Robustly handles cases where bcftools isec finds no shared variants and does not create an output file.

set -euo pipefail

## CONFIGURATION ##
WORKSPACE_DIR="/data/training2/analisis_TFM_Bulida_Precoz/09_samples_double_check"
INPUT_VCF="${WORKSPACE_DIR}/final_bcftools_propDP_15pct.filtered.vcf.gz"
OUTPUT_DIR="${WORKSPACE_DIR}/pairwise_comparison_all"
# Create a dedicated directory for temporary files
TEMP_VCF_DIR="${OUTPUT_DIR}/temp_vcfs"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TEMP_VCF_DIR}"

SAMPLES_TO_CHECK=("mt1" "mt2" "mt3" "mt4" "wt1" "wt2" "wt3" "wt4")
THREADS=8

## SCRIPT LOGIC ##
echo "--- Starting Full Pairwise Sample Identity Check for all 8 Samples ---"

SUMMARY_FILE="${OUTPUT_DIR}/full_pairwise_comparison_summary.txt"
echo -e "Sample1\tSample2\tUnique_To_Sample1\tUnique_To_Sample2\tShared_Variants" > "${SUMMARY_FILE}"

for (( i=0; i<${#SAMPLES_TO_CHECK[@]}; i++ )); do
    for (( j=i+1; j<${#SAMPLES_TO_CHECK[@]}; j++ )); do
        sample1=${SAMPLES_TO_CHECK[i]}
        sample2=${SAMPLES_TO_CHECK[j]}

        echo -e "\n--- Comparing: ${sample1} vs ${sample2} ---"
        
        VCF1="${TEMP_VCF_DIR}/${sample1}.vcf.gz"
        VCF2="${TEMP_VCF_DIR}/${sample2}.vcf.gz"

        # Step 1: Create temporary files for each sample
        bcftools view -s ${sample1} -c1 -Oz -o ${VCF1} ${INPUT_VCF}
        bcftools index ${VCF1}
        bcftools view -s ${sample2} -c1 -Oz -o ${VCF2} ${INPUT_VCF}
        bcftools index ${VCF2}

        # Step 2: Run isec on the temporary files
        PAIR_DIR="${OUTPUT_DIR}/temp_isec_${sample1}_vs_${sample2}"
        # Let bcftools isec generate its standard 4 output directories/files
        bcftools isec -p "${PAIR_DIR}" -c all ${VCF1} ${VCF2}

        # --- THIS IS THE CORRECTED PART ---
        # Step 3: Robustly count the results, assigning 0 if a file doesn't exist
        FILE_UNIQUE1="${PAIR_DIR}/0000.vcf"
        FILE_UNIQUE2="${PAIR_DIR}/0001.vcf"
        FILE_SHARED="${PAIR_DIR}/0002.vcf" # This file contains records from file1 also present in file2

        unique1=$( [ -f "${FILE_UNIQUE1}" ] && grep -vc '^#' "${FILE_UNIQUE1}" || echo 0 )
        unique2=$( [ -f "${FILE_UNIQUE2}" ] && grep -vc '^#' "${FILE_UNIQUE2}" || echo 0 )
        shared=$( [ -f "${FILE_SHARED}" ] && grep -vc '^#' "${FILE_SHARED}" || echo 0 )
        # --- END OF CORRECTION ---

        echo "Results: Unique to ${sample1}: ${unique1} | Unique to ${sample2}: ${unique2} | Shared: ${shared}"
        
        echo -e "${sample1}\t${sample2}\t${unique1}\t${unique2}\t${shared}" >> "${SUMMARY_FILE}"
        
        # Clean up all temporary files for this pair
        rm -r "${PAIR_DIR}"
        rm ${VCF1}* ${VCF2}*
    done
done

# Final cleanup
rmdir "${TEMP_VCF_DIR}"

echo -e "\n### Full pairwise comparison complete! ###"
echo "Summary of results:"
cat "${SUMMARY_FILE}" | column -t