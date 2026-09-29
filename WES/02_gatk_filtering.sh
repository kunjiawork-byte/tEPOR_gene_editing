#!/usr/bin/env bash
# Recovered from supplied WES command notes; analytical options retained.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
gatk GetPileupSummaries \
  -R "$REF" \
  -I "$BASE/mock/mock.aln.bam" \
  -V "$RES/small_exac_common_3.hg38.vcf.gz" \
  -L "$BED" \
  -O "$OUT/mock.pileups.table"
for SAMPLE in 3XBE tEPORCas9; do
gatk LearnReadOrientationModel \
  -I "$OUT/$SAMPLE/$SAMPLE.f1r2.tar.gz" \
  -O "$OUT/$SAMPLE/$SAMPLE.read-orientation-model.tar.gz"
gatk GetPileupSummaries \
  -R "$REF" \
  -I "$BASE/$SAMPLE/$SAMPLE.aln.bam" \
  -V "$RES/small_exac_common_3.hg38.vcf.gz" \
  -L "$BED" \
  -O "$OUT/$SAMPLE/$SAMPLE.pileups.table"
gatk CalculateContamination \
  -I "$OUT/$SAMPLE/$SAMPLE.pileups.table" \
  -matched "$OUT/mock.pileups.table" \
  -O "$OUT/$SAMPLE/$SAMPLE.contamination.table" \
  --tumor-segmentation "$OUT/$SAMPLE/$SAMPLE.segments.table"
gatk FilterMutectCalls \
  -R "$REF" \
  -V "$OUT/$SAMPLE/$SAMPLE.unfiltered.vcf.gz" \
  --contamination-table "$OUT/$SAMPLE/$SAMPLE.contamination.table" \
  --tumor-segmentation "$OUT/$SAMPLE/$SAMPLE.segments.table" \
  --ob-priors "$OUT/$SAMPLE/$SAMPLE.read-orientation-model.tar.gz" \
  -O "$OUT/$SAMPLE/$SAMPLE.filtered.vcf.gz"
done
