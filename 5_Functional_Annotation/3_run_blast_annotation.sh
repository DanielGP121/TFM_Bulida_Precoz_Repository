#!/usr/bin/env bash
# SCRIPT: Local BLASTp of Bulida Precoz proteins vs. Prunus armeniaca reference proteome
# USAGE: ./blastp_prunus_local.sh

set -euo pipefail

#### 1) CONFIGURATION ####
WORKDIR="/data/training2/analisis_TFM_Bulida_Precoz/10_functional_annotation"
INPUT_FASTA="${WORKDIR}/target_proteins_for_analysis.fasta"

# Where we'll put the prunus proteome and BLAST db
REF_FASTA="${WORKDIR}/prunus_armeniaca_refprot.fasta"
DB_PREFIX="${WORKDIR}/prunus_armeniaca_db"

# Output
OUT_FULL="${WORKDIR}/blastp_prunus_armeniaca_local.tsv"
OUT_TOP1="${WORKDIR}/blastp_prunus_armeniaca_local_top1.tsv"

# BLAST parameters
EVALUE="1e-5"
THREADS=8

#### 2) CHECK DEPENDENCIES ####
for tool in wget makeblastdb blastp; do
  if ! command -v ${tool} &> /dev/null; then
    echo "ERROR: ${tool} not found; please install BLAST+ and wget" >&2
    exit 1
  fi
done

#### 3) DOWNLOAD Prunus armeniaca PROTEOME ####
if [ ! -s "${REF_FASTA}" ]; then
  echo "--- Fetching Prunus armeniaca proteome from UniProt ---"
  wget -O "${REF_FASTA}" \
    "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:36596&format=fasta"
else
  echo "--- Reference proteome already present: ${REF_FASTA} ---"
fi

#### 4) MAKE BLAST DATABASE ####
if [ ! -f "${DB_PREFIX}.pin" ]; then
  echo "--- Building BLAST database: ${DB_PREFIX} ---"
  makeblastdb \
    -in "${REF_FASTA}" \
    -dbtype prot \
    -out "${DB_PREFIX}"
else
  echo "--- BLAST DB already exists: ${DB_PREFIX}.* ---"
fi

#### 5) RUN BLASTp LOCALLY ####
echo "--- Running BLASTp (local) of ${INPUT_FASTA} vs ${DB_PREFIX} ---"
blastp \
  -query "${INPUT_FASTA}" \
  -db "${DB_PREFIX}" \
  -evalue "${EVALUE}" \
  -num_threads "${THREADS}" \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -out "${OUT_FULL}"

echo "Full BLASTp results saved to ${OUT_FULL}"

#### 6) EXTRACT TOP‐HIT PER QUERY ####
echo "--- Extracting single best hit per query (highest bitscore) ---"
sort -k1,1 -k6,6nr "${OUT_FULL}" | awk '!seen[$1]++' > "${OUT_TOP1}"
echo "Top‐hit table saved to ${OUT_TOP1}"

#### 7) DONE ####
echo "===== COMPLETE ====="
echo "Full results: ${OUT_FULL}"
echo "Best hits   : ${OUT_TOP1}"


