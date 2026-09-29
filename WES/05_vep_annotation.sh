#!/usr/bin/env bash
# Recovered from supplied WES command notes; analytical options retained.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
# Requires strict VCFs supplied separately: generation command was not recovered.
# The original VEP tabular output options are retained.
for SAMPLE in 3XBE tEPORCas9; do
  test -s "$OUT/$SAMPLE/$SAMPLE.PASS.DP30.ALT5.vcf.gz" || {
    echo "Missing strict VCF for $SAMPLE. See README.md." >&2; exit 1;
  }
done
for SAMPLE in 3XBE tEPORCas9; do

    INPUT="$OUT/$SAMPLE/$SAMPLE.PASS.DP30.ALT5.vcf.gz"
    OUTPUT="$OUT/$SAMPLE/$SAMPLE.VEP.annotation.tsv"

    echo "===== $SAMPLE ====="

    vep \
      --input_file "$INPUT" \
      --output_file "$OUTPUT" \
      --format vcf \
      --tab \
      --cache \
      --offline \
      --dir_cache "$VEP_CACHE" \
      --assembly GRCh38 \
      --species homo_sapiens \
      --fasta "$REF" \
      --symbol \
      --canonical \
      --mane \
      --hgvs \
      --variant_class \
      --pick \
      --fields "Location,Allele,REF_ALLELE,VARIANT_CLASS,SYMBOL,Gene,Feature,BIOTYPE,Consequence,IMPACT,EXON,INTRON,HGVSc,HGVSp,Existing_variation,CANONICAL,MANE_SELECT" \
      --force_overwrite

done
