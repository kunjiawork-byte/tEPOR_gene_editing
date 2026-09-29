#!/usr/bin/env bash
# Recovered from supplied WES command notes; analytical options retained.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
# BED-restricted calling with 100 bp padding, as recorded in uploaded notes.
for SAMPLE in 3XBE tEPORCas9; do
  mkdir -p "$OUT/$SAMPLE"
gatk --java-options "-Xmx16g" Mutect2 \
  -R "$REF" \
  -I "$BASE/$SAMPLE/$SAMPLE.aln.bam" \
  -I "$BASE/mock/mock.aln.bam" \
  -tumor $SAMPLE \
  -normal mock \
  -L "$BED" \
  --interval-padding 100 \
  --germline-resource "$RES/af-only-gnomad.hg38.vcf.gz" \
  --panel-of-normals "$RES/1000g_pon.hg38.vcf.gz" \
  --f1r2-tar-gz "$OUT/$SAMPLE/$SAMPLE.f1r2.tar.gz" \
  -O "$OUT/$SAMPLE/$SAMPLE.unfiltered.vcf.gz" \
  2>&1 | tee "$OUT/$SAMPLE/$SAMPLE_mutect2.log"
done
