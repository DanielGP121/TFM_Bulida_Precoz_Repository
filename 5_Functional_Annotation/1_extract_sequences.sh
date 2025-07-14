#!/usr/bin/env bash
# SCRIPT: Extracts protein sequences for a target list of genes.

set -euo pipefail

## CONFIGURATION ##
WORKSPACE_DIR="/data/training2/analisis_TFM_Bulida_Precoz/10_functional_annotation"
GENOME_FASTA="/data/training2/info/assemblies/BUL_cur_guided.v1.0.fasta"
ANNOTATION_GFF3="/data/training2/info/assemblies/Bulida_annotationRNAseq.gff3"
GENE_LIST_FILE="${WORKSPACE_DIR}/final_gene_list.txt"

# Output files
ALL_PROTEINS_FASTA="${WORKSPACE_DIR}/all_proteins.fasta"
TARGET_PROTEINS_FASTA="${WORKSPACE_DIR}/target_proteins_for_analysis.fasta"

# --- Step 1: Extract ALL protein sequences from the genome ---
if [ ! -s "${ALL_PROTEINS_FASTA}" ]; then # -s checks if file exists AND is not empty
    echo "--- Extracting ALL protein sequences from the genome ---"
    gffread "${ANNOTATION_GFF3}" -g "${GENOME_FASTA}" -y "${ALL_PROTEINS_FASTA}"
else
    echo "--- File with all proteins already exists. Skipping extraction. ---"
fi

# --- Step 2: Select the protein sequences for the genes of interest ---
echo "--- Selecting protein sequences for the target genes ---"
# Use grep to find the headers matching our gene list (-Fwf)
# and print the header line plus the sequence that follows (-A 1)
# This is more robust than seqtk subseq for this task.
grep -A 1 -Fwf "${GENE_LIST_FILE}" "${ALL_PROTEINS_FASTA}" > "${TARGET_PROTEINS_FASTA}"

# Remove the "--" separator lines that grep adds between matches
sed -i '/^--$/d' "${TARGET_PROTEINS_FASTA}"

echo -e "\n### Process complete. ###"
# Verify that the output file now has content
COUNT=$(grep -c ">" "${TARGET_PROTEINS_FASTA}")
echo "Found and extracted ${COUNT} sequences."
echo "The file with the proteins to analyze is: ${TARGET_PROTEINS_FASTA}"