#!/usr/bin/env bash
# tEPOR study: QC for demultiplexed FASTQs; retains FastQC -t 8.
set -euo pipefail
QC_DIR="${QC_DIR:-/path/to/demultiplexed_data}"
cd "$QC_DIR"
# Null-delimited paths allow spaces in filenames and directories.
find . -type f -name 'trimmed-*.fastq.gz' -print0 | xargs -0 fastqc -t 8
multiqc . -o multiqc_report
