#!/usr/bin/env bash
# SCRIPT: Annotate BLASTp hits with UniProt GO terms (accessions endpoint + timeouts)

set -euo pipefail

#### 1) CONFIGURATION ####
WORKDIR="/data/training2/analisis_TFM_Bulida_Precoz/10_functional_annotation"
BLAST_IN="${WORKDIR}/blastp_prunus_armeniaca_local.tsv"
ACCESSIONS_LIST="${WORKDIR}/subject_accessions.list"
GO_MAPPING="${WORKDIR}/uniprot_GO_mapping.tsv"
GO_COLLAPSED="${WORKDIR}/uniprot_GO_mapping.collapsed.tsv"
TMP_BLAST2="${WORKDIR}/blast_for_join.tsv"
FINAL_OUT="${WORKDIR}/blast_with_GO.tsv"

#### 2) EXTRACT UNIQUE ACCESSIONS ####
awk '{ 
  acc=$2; 
  sub(/^tr\|/,"",acc); 
  sub(/\|.*$/,"",acc); 
  print acc 
}' "${BLAST_IN}" | sort -u > "${ACCESSIONS_LIST}"

num_acc=$(wc -l < "${ACCESSIONS_LIST}")
echo "Found ${num_acc} unique UniProt accessions in BLAST hits."

#### 3) FETCH GO IDs VIA ACCESSIONS ENDPOINT ####
accessions=$(paste -sd, "${ACCESSIONS_LIST}")
echo "Querying UniProt accession endpoint for GO_IDs…"

curl -s \
     --connect-timeout 10 \
     -m 30 \
     "https://rest.uniprot.org/uniprotkb/accessions?accessions=${accessions}&fields=accession,go_id&format=tsv" \
  | tail -n +2 > "${GO_MAPPING}" 

if [ ! -s "${GO_MAPPING}" ]; then
  echo "No GO IDs returned. Your proteins may lack UniProt GO cross-references."
  exit 0
fi
echo "Downloaded GO mapping to ${GO_MAPPING}"

#### 4) COLLAPSE MULTIPLE GO_IDs PER ACCESSION ####
awk 'BEGIN{FS=OFS="\t"}{
  acc=$1; gid=$2;
  arr[acc]=(arr[acc]?arr[acc]","gid:gid)
}
END{
  for(a in arr) print a, arr[a]
}' "${GO_MAPPING}" > "${GO_COLLAPSED}"
echo "✔ Collapsed GO_IDs per accession into ${GO_COLLAPSED}"

#### 5) PREPARE BLAST TABLE FOR JOIN ####
awk 'BEGIN{FS=OFS="\t"}{
  acc=$2; sub(/^tr\|/,"",acc); sub(/\|.*$/,"",acc);
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
    $1, acc, $3, $4, $5, $6, $7)
}' "${BLAST_IN}" > "${TMP_BLAST2}"
echo "✔ Prepared BLAST join file ${TMP_BLAST2}"

#### 6) JOIN GO_IDs INTO BLAST HITS ####
echo -e "qseqid\taccession\tpident\tlength\tevalue\tbitscore\ttitle\tGO_IDs" > "${FINAL_OUT}"
join -t $'\t' -1 2 -2 1 \
  <(sort -k2,2 "${TMP_BLAST2}") \
  <(sort -k1,1 "${GO_COLLAPSED}") \
  >> "${FINAL_OUT}"

echo "Annotation complete! See: ${FINAL_OUT}"



