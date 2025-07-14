#!/usr/bin/env bash
# FINAL SCRIPT (Corrected v2): Runs the full downstream analysis pipeline.
# Includes a step to remove the problematic 'AS_FilterStatus' field before merging.

set -euo pipefail

## ============================================================================================
## CONFIGURATION
## ============================================================================================

BASE_DIR="/data/training2/analisis_TFM_Bulida_Precoz"
EXP_DIR="${BASE_DIR}/08_parameter_change_analysis"
FINAL_CANDIDATE_DIR="${EXP_DIR}/final_candidate_sets"
mkdir -p "${FINAL_CANDIDATE_DIR}"

# --- Input Directories ---
INPUT_DIR_DEFAULT="${EXP_DIR}/mutect2_final_workflow/01_mutect2_filtered"
INPUT_DIR_PROP_DP="${EXP_DIR}/mutect2_final_workflow/02_dp_filtered"

# --- Output Prefixes ---
OUTPUT_PREFIX_DEFAULT="mutect2_default_filter"
OUTPUT_PREFIX_PROP_DP="mutect2_plus_propDP_15pct"

# --- General Parameters ---
MUTANT_SAMPLES=("mt1" "mt2" "mt3" "mt4")
THREADS=8
SNPEFF_GENOME_NAME="bulida_v1"
STRICT_FILTER_EXP='N_PASS(GT="alt")==4'


## ============================================================================================
## ANALYSIS FUNCTION
## ============================================================================================

run_analysis_pipeline() {
    local input_dir="$1"
    local output_prefix="$2"
    local temp_dir="${FINAL_CANDIDATE_DIR}/temp_${output_prefix}"
    mkdir -p "${temp_dir}"

    echo "======================================================="
    echo "--- Starting Downstream Analysis for prefix: ${output_prefix} ---"
    echo "--- Input directory: ${input_dir} ---"

    # --- Step 1: Pre-processing before merge ---
    local VCF_LIST_FINAL=()
    for sample in "${MUTANT_SAMPLES[@]}"; do
        local original_vcf=$(find "${input_dir}" -name "${sample}.*.vcf.gz" | head -n 1)
        local single_sample_vcf="${temp_dir}/${sample}.single.vcf.gz"
        local clean_vcf="${temp_dir}/${sample}.clean.vcf.gz"
        
        if [ -f "$original_vcf" ]; then
            echo "Processing ${sample}:"
            # 1a: Extract only the mutant sample to solve duplicate name issue
            echo "  - Extracting mutant sample..."
            bcftools view --threads "${THREADS}" --samples "${sample}" -Oz -o "${single_sample_vcf}" "${original_vcf}"
            
            # 1b: Remove the problematic INFO field to solve merge error
            echo "  - Removing AS_FilterStatus field..."
            bcftools annotate --threads "${THREADS}" -x INFO/AS_FilterStatus -Oz -o "${clean_vcf}" "${single_sample_vcf}"
            
            bcftools index --threads "${THREADS}" "${clean_vcf}"
            VCF_LIST_FINAL+=("${clean_vcf}")
        else
            echo "ERROR: Could not find VCF for ${sample} in ${input_dir}"
            exit 1
        fi
    done

    # 2. Merge the super-clean, single-sample VCFs
    echo "Step 2: Merging 4 clean single-sample VCFs..."
    local MERGED_VCF="${temp_dir}/merged.vcf.gz"
    bcftools merge --threads "${THREADS}" -Oz -o "${MERGED_VCF}" "${VCF_LIST_FINAL[@]}"
    bcftools index --threads "${THREADS}" "${MERGED_VCF}"

    # 3. Apply the strict biological filter (4 out of 4)
    local STRICT_CANDIDATES_VCF="${FINAL_CANDIDATE_DIR}/${output_prefix}_strict_candidates.vcf.gz"
    echo "Step 3: Applying strict biological filter (${STRICT_FILTER_EXP})..."
    bcftools view --threads "${THREADS}" --include "${STRICT_FILTER_EXP}" -Oz -o "${STRICT_CANDIDATES_VCF}" "${MERGED_VCF}"
    bcftools index --threads "${THREADS}" "${STRICT_CANDIDATES_VCF}"
    local COUNT=$(bcftools view -H "${STRICT_CANDIDATES_VCF}" | wc -l)
    echo "Found ${COUNT} strict (4/4) candidates for this filter set."

    # 4. Annotate and Extract Results
    if [ "${COUNT}" -gt 0 ]; then
        local ANNOTATED_VCF="${FINAL_CANDIDATE_DIR}/${output_prefix}_final.ann.vcf"
        local TABLE_OUTPUT="${FINAL_CANDIDATE_DIR}/${output_prefix}_final_table.tsv"
        local GENE_LIST_OUTPUT="${FINAL_CANDIDATE_DIR}/${output_prefix}_final_gene_list.txt"

        echo "Step 4: Annotating and extracting final results..."
        snpEff ann -v "${SNPEFF_GENOME_NAME}" "${STRICT_CANDIDATES_VCF}" > "${ANNOTATED_VCF}"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' "${ANNOTATED_VCF}" > "${TABLE_OUTPUT}"
        cut -f 5 "${TABLE_OUTPUT}" | cut -d '|' -f 5 | sort -u | grep -v '^\s*$' > "${GENE_LIST_OUTPUT}"
    else
        echo "No candidates found, skipping annotation and extraction."
    fi
    
    # 5. Cleanup
    rm -rf "${temp_dir}"
    echo "Analysis for prefix '${output_prefix}' is complete."
}


## ============================================================================================
## SCRIPT EXECUTION
## ============================================================================================

# Execute the pipeline for the first filter strategy
run_analysis_pipeline "${INPUT_DIR_DEFAULT}" "${OUTPUT_PREFIX_DEFAULT}"

# Execute the pipeline for the second filter strategy
run_analysis_pipeline "${INPUT_DIR_PROP_DP}" "${OUTPUT_PREFIX_PROP_DP}"

echo -e "\n### Both analysis pipelines have finished successfully! ###"
echo "Final results are in: ${FINAL_CANDIDATE_DIR}"