#!/usr/bin/env bash
set -euo pipefail

# concat_wt4.sh
# Concatenate the two A_2 lane files (wt4) in 00_raw_data into single R1/R2 FASTQs

WORKDIR="/data/training2/analisis_TFM_Bulida_Precoz/00_raw_data"
OUT_R1="$WORKDIR/wt4_R1.fq.gz"
OUT_R2="$WORKDIR/wt4_R2.fq.gz"

# Debug: show what will be run
echo "=== DEBUG ==="
echo "Concatenating lane-1: $WORKDIR/A_2*_L2_1.fq.gz → $OUT_R1"
echo "Concatenating lane-2: $WORKDIR/A_2*_L2_2.fq.gz → $OUT_R2"
echo "============="

echo "Starting concatenation of wt4 (A_2) reads…"

# Lane-1 → R1
echo "zcat $WORKDIR/A_2*_L2_1.fq.gz | gzip > $OUT_R1"
zcat "$WORKDIR"/A_2*_L2_1.fq.gz | gzip > "$OUT_R1"

# Lane-2 → R2
echo "zcat $WORKDIR/A_2*_L2_2.fq.gz | gzip > $OUT_R2"
zcat "$WORKDIR"/A_2*_L2_2.fq.gz | gzip > "$OUT_R2"

echo "==== Concatenation of wt4 reads complete ===="
