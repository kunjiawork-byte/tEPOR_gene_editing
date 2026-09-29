#!/usr/bin/env bash
# Recovered from supplied WES command notes; analytical options retained.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
# Strict set only: PASS, tumor DP >= 30, first ALT read count >= 5.
# Original AD[1]/AD[2] handling and SNV/INDEL classification are retained.
# Multiallelic and MNV handling require the missing later normalization stage.
for SAMPLE in 3XBE tEPORCas9; do

    VCF="$OUT/$SAMPLE/$SAMPLE.filtered.vcf.gz"
    TSV="$OUT/$SAMPLE/$SAMPLE.PASS.DP30.ALT5.tsv"

    echo "===== $SAMPLE ====="

    bcftools view -s "$SAMPLE" -f PASS "$VCF" -Ou |
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%DP\t%AD]\n' |
    awk -v OFS="\t" '
    BEGIN {
        print "CHROM","POS","REF","ALT","FILTER","DP",
              "REF_READS","ALT_READS","VAF","TYPE"
    }
    {
        split($7, ad, ",")
        refreads = ad[1]
        altreads = ad[2]

        if ($6 >= 30 && altreads >= 5) {

            vaf = altreads / (refreads + altreads)

            if (length($3)==1 && length($4)==1)
                type="SNV"
            else
                type="INDEL"

            print $1,$2,$3,$4,$5,$6,
                  refreads,altreads,vaf,type
        }
    }' > "$TSV"

done
