#!/usr/bin/env bash
set -euo pipefail
trap 'echo "ERROR: Script failed on line $LINENO with exit code $?"; exit 1' ERR

# run_haplotypecaller_filtering.sh
# This script performs joint-calling on GVCFs from HaplotypeCaller,
# applies technical quality filters, and then performs biological filtering
# to identify candidate somatic mutations.

# ── CONFIGURATION ────────────────────────────────────────────────────────────────
# Directory containing the input GVCF files
GVCF_DIR="/data/training2/analisis_TFM_Bulida_Precoz/04_variant_calling_gatk"

# Main output directory for this filtering workflow
FILTER_DIR="${GVCF_DIR}/filtering_workflow"

# Reference Genome FASTA file
REF_GENOME="/data/training2/info/assemblies/BUL_cur_guided.v1.0.fasta"

# GATK and Java options
MEM_MB=20480 # Approx 20GB
THREADS=8    # Number of threads for parallel processes

# Sample lists
WILD_TYPE_SAMPLES="wt1,wt2,wt3,wt4"
MUTANT_SAMPLES="mt1,mt2,mt3,mt4"

# ── PREPARATION ───────────────────────────────────────────────────────────────────
echo "INFO: GATK Joint-Calling and Filtering script started at $(date)"
mkdir -p "$FILTER_DIR" "${FILTER_DIR}/genomicsdb" "${FILTER_DIR}/logs"

# --- STAGE 1: JOINT-CALLING ───────────────────────────────────────────────────────
# Consolidate GVCFs into a GenomicsDB workspace for efficient joint-calling.
DB_PATH="${FILTER_DIR}/genomicsdb/bulida_cohort_db"
RAW_VCF="${FILTER_DIR}/bulida_cohort.raw.vcf.gz"
INTERVAL_LIST="${FILTER_DIR}/contigs.list"

echo "--- Stage 1: Joint-Calling ---"
if [ ! -d "$DB_PATH" ]; then
    echo "  Step 1/2: Consolidating GVCFs with GenomicsDBImport..."
    GVCF_ARGS=()
    for SAMPLE in $(echo "$WILD_TYPE_SAMPLES,$MUTANT_SAMPLES" | tr ',' ' '); do
        GVCF_ARGS+=(--variant "${GVCF_DIR}/${SAMPLE}.g.vcf.gz")
    done
    
    # Create interval list from reference genome index
    cut -f1 "${REF_GENOME}.fai" > "$INTERVAL_LIST"
    echo "  INFO: Created interval list at ${INTERVAL_LIST}"

    gatk --java-options "-Xmx${MEM_MB}m" GenomicsDBImport \
        "${GVCF_ARGS[@]}" \
        --genomicsdb-workspace-path "$DB_PATH" \
        -L "$INTERVAL_LIST" \
        --reader-threads ${THREADS} \
        --genomicsdb-shared-posixfs-optimizations
else
    echo "  INFO: GenomicsDB workspace already exists. Skipping import."
fi

# Perform joint genotyping
if [ ! -f "$RAW_VCF" ]; then
    echo "  Step 2/2: Joint-genotyping with GenotypeGVCFs..."
    gatk --java-options "-Xmx${MEM_MB}m" GenotypeGVCFs \
        -R "$REF_GENOME" \
        -V "genomicsdb://${DB_PATH}" \
        -O "$RAW_VCF"
else
    echo "  INFO: Raw joint-called VCF already exists. Skipping GenotypeGVCFs."
fi

# --- STAGE 2: TECHNICAL QUALITY FILTERING ───────────────────────────────────────
echo -e "\n--- Stage 2: Technical Quality Filtering ---"
FILTERED_SNPS_VCF="${FILTER_DIR}/bulida_cohort.filtered_snps.vcf.gz"
FILTERED_INDELS_VCF="${FILTER_DIR}/bulida_cohort.filtered_indels.vcf.gz"

if [ ! -f "$FILTERED_SNPS_VCF" ] || [ ! -f "$FILTERED_INDELS_VCF" ]; then
    # Step 2.1: Select SNPs and Indels into separate files
    RAW_SNPS_VCF="${FILTER_DIR}/bulida_cohort.raw_snps.vcf.gz"
    RAW_INDELS_VCF="${FILTER_DIR}/bulida_cohort.raw_indels.vcf.gz"
    echo "  Step 1/4: Selecting SNPs..."
    gatk SelectVariants -V "$RAW_VCF" -select-type SNP -O "$RAW_SNPS_VCF"
    echo "  Step 2/4: Selecting Indels..."
    gatk SelectVariants -V "$RAW_VCF" -select-type INDEL -O "$RAW_INDELS_VCF"

    # Step 2.2: Apply hard filters
    echo "  Step 3/4: Applying filters to SNPs..."
    gatk VariantFiltration \
        -V "$RAW_SNPS_VCF" \
        -O "$FILTERED_SNPS_VCF" \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "SNP_HARD_FILTER"

    echo "  Step 4/4: Applying filters to Indels..."
    gatk VariantFiltration \
        -V "$RAW_INDELS_VCF" \
        -O "$FILTERED_INDELS_VCF" \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "INDEL_HARD_FILTER"
    
    echo "  INFO: Technical filtering complete. Filtered SNPs and Indels are in separate files."

else
    echo "  INFO: Filtered VCFs already exist. Skipping technical filtering."
fi

# --- STAGE 3: BIOLOGICAL FILTERING FOR CANDIDATE MUTATIONS ───────────────────
echo -e "\n--- Stage 3: Biological Filtering for Candidate Mutations ---"
CANDIDATES_VCF="${FILTER_DIR}/candidate_somatic_mutations.vcf"

# Construct bcftools filter expression
# All WT must be homozygous reference (0/0)
# All MT must NOT be homozygous reference (^0/0)
WT_EXPR=$(echo "$WILD_TYPE_SAMPLES" | tr ',' ' ' | sed 's/\([^ ]*\)/GT[\"&"]=\"0\/0\"/g' | tr ' ' '&')
MT_EXPR=$(echo "$MUTANT_SAMPLES" | tr ',' ' ' | sed 's/\([^ ]*\)/GT[\"&"]!=\"0\/0\"/g' | tr ' ' '&')
FULL_EXPR="${WT_EXPR} & ${MT_EXPR}"

echo "  Filtering with expression: ${FULL_EXPR}"

CANDIDATE_SNPS="${FILTER_DIR}/candidate_snps.vcf.gz"
CANDIDATE_INDELS="${FILTER_DIR}/candidate_indels.vcf.gz"

# Filter SNPs that passed technical filters AND meet biological criteria
bcftools view -i 'FILTER="PASS"' "$FILTERED_SNPS_VCF" | bcftools view -i "$FULL_EXPR" -Oz -o "$CANDIDATE_SNPS"

# Filter Indels that passed technical filters AND meet biological criteria
bcftools view -i 'FILTER="PASS"' "$FILTERED_INDELS_VCF" | bcftools view -i "$FULL_EXPR" -Oz -o "$CANDIDATE_INDELS"

# Combine candidate SNPs and Indels into the final file
tabix -p vcf "$CANDIDATE_SNPS"
tabix -p vcf "$CANDIDATE_INDELS"
bcftools concat -a "$CANDIDATE_SNPS" "$CANDIDATE_INDELS" -o "$CANDIDATES_VCF"

echo "============================================================"
echo "Biological filtering complete."
echo "Final candidate mutations are in: $CANDIDATES_VCF"
echo "Script finished at $(date)"
