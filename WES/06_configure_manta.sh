#!/usr/bin/env bash
# Configuration command recovered from the header of uploaded runWorkflow.py.
# This configures the recorded 3XBE single-sample exome workflow; it does not run it.
set -euo pipefail
MANTA_BAM="${MANTA_BAM:-/path/to/3XBE.manta.clean.bam}"
MANTA_REF="${MANTA_REF:-/path/to/GRCh38.canonical.fa}"
MANTA_RUN_DIR="${MANTA_RUN_DIR:-/path/to/manta/3XBE}"
configManta.py --bam "$MANTA_BAM" --referenceFasta "$MANTA_REF" --runDir "$MANTA_RUN_DIR" --exome
# Use the generated runWorkflow.py in its own run directory with Python 2.
# Original run flags and downstream deletion classification were not supplied.
