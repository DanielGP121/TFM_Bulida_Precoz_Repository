#!/usr/bin/env bash
# SCRIPT: Creates a comprehensive gene -> GO mapping file for the entire proteome.

set -euo pipefail

#### CONFIGURATION ####
WORKDIR="/data/training2/analisis_TFM_Bulida_Precoz/10_functional_annotation"
ALL_PROTEINS_FASTA="${WORKDIR}/all_proteins.fasta"
ALL_PROTEINS_CLEANED_FASTA="${WORKDIR}/all_proteins.cleaned.fasta"
DB_PREFIX="${WORKDIR}/prunus_armeniaca_db"
UNIVERSE_FILE_OUT="${WORKDIR}/gene2go_universe_full.tsv"

# Temporary files
BLAST_OUT_FULL="${WORKDIR}/temp_blast_all_proteins.tsv"
ACCESSIONS_LIST_FULL="${WORKDIR}/temp_all_accessions.list"
GO_MAPPING_FULL="${WORKDIR}/temp_all_uniprot_go_mapping.tsv"

THREADS=8
BATCH_SIZE=400 # Number of accessions to query per request

#### SCRIPT LOGIC ####
echo "--- Creating Full GO Universe for Enrichment Analysis ---"

# --- Step 0: Clean the protein FASTA file ---
if [ ! -s "${ALL_PROTEINS_CLEANED_FASTA}" ]; then
    echo "Step 0: Cleaning the protein FASTA file..."
    awk '/^>/ {print; next} {gsub(/\*|X/,""); print}' "${ALL_PROTEINS_FASTA}" > "${ALL_PROTEINS_CLEANED_FASTA}"
    echo "Cleaned FASTA created."
else
    echo "Step 0: Cleaned FASTA file already exists. Skipping."
fi

# --- Step 1: Run BLASTp ---
if [ ! -s "${BLAST_OUT_FULL}" ]; then
    echo "Step 1: Running BLASTp for all proteins..."
    blastp \
        -query "${ALL_PROTEINS_CLEANED_FASTA}" \
        -db "${DB_PREFIX}" \
        -evalue 1e-5 \
        -num_threads "${THREADS}" \
        -outfmt "6 qseqid sseqid" \
        -out "${BLAST_OUT_FULL}"
else
    echo "Step 1: BLASTp output already exists. Skipping."
fi

# --- Step 2: Extract UniProt accessions ---
if [ ! -s "${ACCESSIONS_LIST_FULL}" ]; then
    echo "Step 2: Extracting UniProt accessions..."
    awk '{ acc=$2; sub(/^tr\|/,"",acc); sub(/\|.*$/,"",acc); print acc }' "${BLAST_OUT_FULL}" \
    | sort -u > "${ACCESSIONS_LIST_FULL}"
else
    echo "Step 2: Accession list already exists. Skipping."
fi

# --- Step 3: Fetch GO terms from UniProt in batches (THE FIX) ---
if [ ! -s "${GO_MAPPING_FULL}" ]; then
    echo "Step 3: Fetching GO terms from UniProt in batches..."
    # Create the output file so we can append to it
    touch "${GO_MAPPING_FULL}"
    
    # Read the accession list file line by line and process in batches
    mapfile -t accessions_array < "${ACCESSIONS_LIST_FULL}"
    total_accessions=${#accessions_array[@]}
    
    for (( i=0; i<total_accessions; i+=BATCH_SIZE )); do
        # Get a batch of accessions
        batch_array=("${accessions_array[@]:i:BATCH_SIZE}")
        batch_string=$(IFS=,; echo "${batch_array[*]}")
        
        echo "  - Fetching batch starting at accession #${i}..."
        
        # Run curl for the batch and append (>>) to the output file
        curl -s --connect-timeout 20 -m 60 \
             "https://rest.uniprot.org/uniprotkb/accessions?accessions=${batch_string}&fields=accession,go_id&format=tsv" \
        | tail -n +2 >> "${GO_MAPPING_FULL}"

        sleep 1 # Be nice to the UniProt server and wait 1 second
    done
else
    echo "Step 3: GO mapping file already exists. Skipping."
fi

# --- Step 4: Create the final two-column gene2go file ---
if [ ! -s "${UNIVERSE_FILE_OUT}" ]; then
    echo "Step 4: Creating the final universe file..."
    awk 'BEGIN{FS=OFS="\t"} FNR==NR{go[$1]=$2; next} { acc=$2; sub(/^tr\|/,"",acc); sub(/\|.*$/,"",acc); if(acc in go) print $1, go[acc] }' \
        <(awk 'BEGIN{FS=OFS="\t"}{split($2, go_list, "; "); for(i in go_list) print $1, go_list[i]}' "${GO_MAPPING_FULL}") \
        "${BLAST_OUT_FULL}" \
    | awk 'BEGIN{FS=OFS="\t"}{split($2, go_list, ","); for(i in go_list) print $1, go_list[i]}' \
    | sort -u > "${UNIVERSE_FILE_OUT}"
else
    echo "Step 4: Final universe file already exists. Skipping."
fi

echo -e "\n### Universe file creation workflow finished. ###"
echo "Check the final output at: ${UNIVERSE_FILE_OUT}"