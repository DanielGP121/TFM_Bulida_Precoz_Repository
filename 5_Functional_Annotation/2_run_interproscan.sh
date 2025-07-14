#!/usr/bin/env bash
# SCRIPT: Runs InterProScan on the target protein sequences.

set -euo pipefail

# 1) Determine Conda and InterProScan installation paths
PREFIX="${CONDA_PREFIX}"
IPS_HOME="${PREFIX}/share/InterProScan"

# 2) Index HMMER libraries if the .h3* files are missing
for lib in \
  "data/superfamily/1.75/hmmlib_1.75" \
  "data/pirsf/3.10/sf_hmm_all"; do

  full_path="${IPS_HOME}/${lib}"
  if [[ -f "${full_path}" && ! -f "${full_path}.h3f" ]]; then
    echo "Indexing ${full_path} with hmmpress..."
    "${PREFIX}/bin/hmmpress" "${full_path}"
  fi
done

## CONFIGURATION ##
WORKSPACE_DIR="/data/training2/analisis_TFM_Bulida_Precoz/10_functional_annotation"
INPUT_FASTA="${WORKSPACE_DIR}/target_proteins_for_analysis.fasta"
OUTPUT_FILE_PREFIX="${WORKSPACE_DIR}/interproscan_results"

# Number of CPUs to use
CPU=8

echo "--- Starting InterProScan analysis ---"
echo "Input FASTA: ${INPUT_FASTA}"
echo "This may take several hours..."

# Run InterProScan with TSV output, GO terms, and pathway annotation
interproscan.sh \
    -i "${INPUT_FASTA}" \
    -f tsv \
    -o "${OUTPUT_FILE_PREFIX}.tsv" \
    --iprlookup \
    --goterms \
    --pathways \
    --disable-precalc \
    -cpu "${CPU}"

echo -e "\n### InterProScan analysis completed. ###"
echo "Results written to: ${OUTPUT_FILE_PREFIX}.tsv"

