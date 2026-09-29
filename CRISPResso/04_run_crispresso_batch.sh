#!/usr/bin/env bash
# tEPOR study: retain the four CRISPRessoBatch parameter sets in the source notes.
set -euo pipefail
BATCH_SETTINGS="${BATCH_SETTINGS:-CRISPResso_batch_input.txt}"
MODE="${1:?Usage: bash 04_run_crispresso_batch.sh cbe|cbe_ga|cas9_cut|cas9_wide}"
options=()
case "$MODE" in
  cbe) options=(--quantification_window_center -10 --quantification_window_size 20 --base_editor_output) ;;
  cbe_ga) options=(--quantification_window_center -10 --quantification_window_size 20 --base_editor_output --conversion_nuc_from G --conversion_nuc_to A) ;;
  cas9_cut) options=(--quantification_window_center -3 --quantification_window_size 1) ;;
  cas9_wide) options=(--quantification_window_center -10 --quantification_window_size 20) ;;
  *) echo "Unknown mode: $MODE" >&2; exit 1 ;;
esac
CRISPRessoBatch --batch_settings "$BATCH_SETTINGS" "${options[@]}" --n_processes 4
