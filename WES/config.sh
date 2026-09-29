#!/usr/bin/env bash
# tEPOR gene editing study: shared WES paths.
# Source this file, or export these variables before running each stage.
BASE="${BASE:-/path/to/WESsample/samples}"
REF="${REF:-/path/to/GRCh38.primary_assembly.genome.fa}"
RES="${RES:-/path/to/gatk_somatic_hg38}"
BED="${BED:-/path/to/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed}"
OUT="${OUT:-/path/to/WESsample/mutect2_reanalysis}"
VEP_CACHE="${VEP_CACHE:-/path/to/vep_cache}"
export BASE REF RES BED OUT VEP_CACHE
