# WES analysis: recovered code subset

This folder contains commands recovered from the supplied WES notes and two Python files. This is a partial reconstruction, not the complete final manuscript workflow. Resource locations and software environments must be configured before use. Analytical parameters were retained; repeated per-sample commands were consolidated into loops.

## Files and execution

Activate an environment providing GATK, bcftools, and Python 3. Edit `config.sh` or export its variables. Run from any directory:

```bash
bash 01_mutect2.sh
bash 02_gatk_filtering.sh
bash 03_strict_variant_tables.sh
source config.sh
python3 04_analyze_snv_spectrum.py
```

`01_mutect2.sh` starts with vendor-aligned BAMs for 3XBE, tEPORCas9, and matched normal mock. It retains the supplied Twist target BED and 100 bp padding. The reference FASTA, BAM contigs, resource contigs, and sample names must match. The earlier chromosome-restricted retry was mentioned in previous context but not included in these uploaded notes; no chromosome restriction was invented here.

`02_gatk_filtering.sh` includes mock and tumor pileups, orientation models, matched-normal contamination estimates, segmentation, and FilterMutectCalls.

`03_strict_variant_tables.sh` reproduces the original first-ALT handling and emits `Sample.PASS.DP30.ALT5.tsv`. It does not normalize or split multiallelic variants. TYPE labels follow the source: single-base REF and ALT are SNV; everything else is INDEL. VAF is ALT_READS / (REF_READS + ALT_READS), not FORMAT/AF.

`04_analyze_snv_spectrum.py` reads those strict tables and writes variant classification, raw 12-class spectrum, pyrimidine-centered 6-class spectrum, and C>T/G>A summary. It assumes nonempty SNV sets as in the original source. These are strict-set outputs, not the later FINAL 51-variant summaries.

`05_vep_annotation.sh` preserves the recovered offline GRCh38 VEP tabular annotation command. It requires pre-existing strict VCFs. Their generation command was not supplied, so this stage is not directly connected to stage 03. VEP and cache versions were not recorded in the supplied files.

`06_configure_manta.sh` retains only the 3XBE single-BAM `--exome` configuration recorded in the auto-generated Manta launcher. Configure Manta in its Python 2 environment with the cleaned BAM and canonical reference. Their generation, workflow run flags, completion, and downstream review are not established by the uploaded launcher. Regenerate the launcher locally instead of copying the vendor-generated Python file containing fixed installation and pickle paths.

## Recovered results and missing final stages

The source reports the strict set as 3XBE 32 variants (31 SNV, 1 INDEL) and tEPORCas9 13 variants (10 SNV, 3 INDEL). These numbers are historical QC targets, not hardcoded filters.

The later expanded set reportedly contains 34 + 17 = 51 normalized variants. The complex-event rescue, final normalization, corrected annotation merge, final consequence categories, FINAL_* summary generation, FINAL_DEL_1kb_classification_clean.tsv generation, and A-C plotting code were not supplied. This subset cannot reproduce those final results.

The uploaded legacy VEP merge used (chromosome, VEP start, VEP allele) to match (CHROM, POS, ALT). That key is unreliable for indels because VEP can change coordinates and trim alleles. The notes indicate later correction of affected annotations but do not include that corrected code. The legacy merge was therefore documented rather than included as a recommended executable.

The last consequence table differs from previously reported final category assignments (for example, intronic/noncoding and upstream/downstream categories). Without the corrected merge and final classification code, it cannot establish the final functional summary.

## Software

Tools in recovered commands: GATK, bcftools, Python 3, offline Ensembl VEP, and Manta with Python 2. Record the actual versions and reference/cache provenance used for the study before final publication.

## Validation

Shell syntax and Python syntax can be checked locally. Full execution requires the original BAMs, reference/resources, VEP cache, and tools; no full data rerun was performed while preparing this folder.
