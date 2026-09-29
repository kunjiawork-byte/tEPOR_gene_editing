#!/usr/bin/env bash
# tEPOR study: paired-end demultiplexing using anchored R2 barcodes.
# Retains the deepseq block's -e 0.1 and -q 20,20.
set -euo pipefail
BASE_DIR="${BASE_DIR:-/path/to/00_fastq}"
SAMPLE_PATTERN="${SAMPLE_PATTERN:-KJ28-D4-*}"
shopt -s nullglob
found=0
for SAMPLE_DIR in "$BASE_DIR"/$SAMPLE_PATTERN; do
  [[ -d "$SAMPLE_DIR" ]] || continue
  found=1
  r1_files=("$SAMPLE_DIR"/*_R1_001.fastq.gz)
  r2_files=("$SAMPLE_DIR"/*_R2_001.fastq.gz)
  if (( ${#r1_files[@]} != 1 || ${#r2_files[@]} != 1 )); then
    echo "Expected exactly one R1/R2 file in $SAMPLE_DIR" >&2; exit 1
  fi
  BARCODE_FILE="$SAMPLE_DIR/barcodes.fasta"
  [[ -s "$BARCODE_FILE" ]] || { echo "Missing $BARCODE_FILE" >&2; exit 1; }
  OUTPUT_DIR="$SAMPLE_DIR/demux"
  mkdir -p "$OUTPUT_DIR"
  cutadapt -e 0.1 -q 20,20 -g ^file:"$BARCODE_FILE" \
    -o "$OUTPUT_DIR/trimmed-{name}.R2.fastq.gz" \
    -p "$OUTPUT_DIR/trimmed-{name}.R1.fastq.gz" \
    "${r2_files[0]}" "${r1_files[0]}" > "$OUTPUT_DIR/cutadapt.log"
done
(( found )) || { echo "No sample directories matched $SAMPLE_PATTERN" >&2; exit 1; }
