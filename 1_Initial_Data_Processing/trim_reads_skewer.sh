#!/usr/bin/env bash
set -euo pipefail
trap 'echo "Error on line $LINENO"; exit 1' ERR

# trim_reads_skewer.sh
# Trim paired‐end FASTQ files using skewer for all samples wt1…mt4

# ── CONFIGURATION ────────────────────────────────────────────────────────────────
RAWDIR="/data/training2/analisis_TFM_Bulida_Precoz/00_raw_data"
TRIMDIR="/data/training2/analisis_TFM_Bulida_Precoz/01_trimming"
ADAPTERS="/data/training2/info/shqMergedAdapters_Primers_representative_rc.fa"
THREADS=8

# ── CHECK PREREQUISITES ──────────────────────────────────────────────────────────
command -v skewer >/dev/null || { echo "ERROR: skewer not found in PATH"; exit 1; }
[[ -d "$RAWDIR" ]]    || { echo "ERROR: RAWDIR '$RAWDIR' does not exist"; exit 1; }
mkdir -p "$TRIMDIR"
[[ -w "$TRIMDIR" ]]   || { echo "ERROR: cannot write to TRIMDIR '$TRIMDIR'"; exit 1; }
[[ -r "$ADAPTERS" ]]  || { echo "ERROR: adapter file '$ADAPTERS' not readable"; exit 1; }

# ── DEBUG INFO ──────────────────────────────────────────────────────────────────
echo "===== DEBUG: verifying paths & settings ====="
echo "Script   : $0"
echo "RAWDIR   : $RAWDIR"
echo "TRIMDIR  : $TRIMDIR"
echo "ADAPTERS : $ADAPTERS"
echo "THREADS  : $THREADS"
echo "============================================="

# ── TRIMMING LOOP ───────────────────────────────────────────────────────────────
for SAMPLE in wt1 wt2 wt3 wt4 mt1 mt2 mt3 mt4; do
  echo "Trimming sample $SAMPLE"

  # Gather all matching R1/R2 files
  R1_FILES=( "$RAWDIR"/"${SAMPLE}"_R1*.gz )
  R2_FILES=( "$RAWDIR"/"${SAMPLE}"_R2*.gz )

  # Skip if none found
  if [[ ! -e "${R1_FILES[0]}" ]]; then
    echo "No R1 files for $SAMPLE, skipping."
    continue
  fi
  if [[ ! -e "${R2_FILES[0]}" ]]; then
    echo "No R2 files for $SAMPLE, skipping."
    continue
  fi

  echo " Found R1: ${R1_FILES[*]}"
  echo " Found R2: ${R2_FILES[*]}"

  # Run skewer
  if skewer \
      -r 0.1 \
      -d 0.05 \
      -k 8 \
      -q 20 \
      -l 75 \
      -m pe \
      -t "$THREADS" \
      -x "$ADAPTERS" \
      "${R1_FILES[@]}" \
      "${R2_FILES[@]}" \
      -z \
      -o "$TRIMDIR/$SAMPLE"
  then
    echo "$SAMPLE trimmed successfully"
  else
    echo "Skewer failed on $SAMPLE" >&2
  fi
done

echo "All samples processed."
