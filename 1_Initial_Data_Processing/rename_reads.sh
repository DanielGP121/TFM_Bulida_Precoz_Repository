#!/usr/bin/env bash
set -euo pipefail
# rename_reads.sh
# Rename FASTQ reads to wt1...mt4 
# skipping A_2 (wt4) reads which were already renamed.

WORKDIR="/data/training2/analisis_TFM_Bulida_Precoz/00_raw_data"

# Mapping: original filename → new filename
declare -A MAP=(
  [E_R1_combined.fastq.gz]=wt1_R1.fastq.gz
  [E_R2_combined.fastq.gz]=wt1_R2.fastq.gz
  [F_R1_combined.fastq.gz]=wt2_R1.fastq.gz
  [F_R2_combined.fastq.gz]=wt2_R2.fastq.gz
  [G_R1_combined.fastq.gz]=mt1_R1.fastq.gz
  [G_R2_combined.fastq.gz]=mt1_R2.fastq.gz
  [H_R1_combined.fastq.gz]=mt2_R1.fastq.gz
  [H_R2_combined.fastq.gz]=mt2_R2.fastq.gz

  [A_1_BDPL200002386-1A_HKHF3DSXY_L2_1.fq.gz]=wt3_R1.fq.gz
  [A_1_BDPL200002386-1A_HKHF3DSXY_L2_2.fq.gz]=wt3_R2.fq.gz
  # A_2 reads have been renamed already, so they're omitted here.
  [A_3_BDPL200002388-1A_HKHF3DSXY_L2_1.fq.gz]=mt3_R1.fq.gz
  [A_3_BDPL200002388-1A_HKHF3DSXY_L2_2.fq.gz]=mt3_R2.fq.gz
  [A_4_BDPL200002389-1A_HKHF3DSXY_L2_1.fq.gz]=mt4_R1.fq.gz
  [A_4_BDPL200002389-1A_HKHF3DSXY_L2_2.fq.gz]=mt4_R2.fq.gz
)

echo "==== Renaming reads in $WORKDIR ===="
for orig in "${!MAP[@]}"; do
  src="$WORKDIR/$orig"
  dst="$WORKDIR/${MAP[$orig]}"
  if [[ -e "$src" ]]; then
    echo "Renaming $orig → ${MAP[$orig]}"
    mv --verbose "$src" "$dst"
  else
    echo "Warning: $orig not found, skipping."
  fi
done
echo "==== Renaming complete ===="
