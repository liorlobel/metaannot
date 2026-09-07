#!/usr/bin/env bash
# The runs. CLAUDE.md rule 6: do this inside tmux — an SSH drop mid-run leaves partial
# state. It resumes correctly, but only if the process was allowed to record what it
# finished.
#     tmux new -s metaannot
#     ./02_run.sh 2>&1 | tee run_$(date +%F_%H%M).log
#
# Sequential on purpose. Each dataset already uses 22 cpu / 80 GB across 3 concurrent
# stages; running datasets in parallel on top of that would oversubscribe the box, and
# four hmmsearch jobs against Pfam-A alongside InterProScan is the usual squeeze.
#
# If a dataset dies, rerun this script: finished stages are cached and skipped, only the
# failed one and its dependents rerun. Do NOT add a bare --force.
set -uo pipefail
cd "$(dirname "$0")"
MA=${MA:-$HOME/metaannot/metaannot.py}
PY=${PY:-python3}
ORDER="DDA_mouse_microbiome_nitrosylation DDA_mouse_microbiome_acetylation mouse_dietary_intervnetion Heyer_IBD_metaproteomics DDA_mouse_microbiome_acetylation2 UC_nitsan_2023_1st_batch Acute_colitis_cohort Results_Combined_UC_metaprotoemics"
for d in $ORDER; do
  [ -f "$d/config.yaml" ] || { echo "!! no config for $d"; continue; }
  echo; echo "################ $d  ($(date +%H:%M)) ################"
  "$PY" "$MA" run --config "$d/config.yaml" 2>&1 | tail -80
  echo "---- $d bins ----"; column -t "$d/results/bin_summary.tsv" 2>/dev/null || echo "(no bin_summary — check the log)"
done
echo; echo "grep -h 'invisible to KEGG' */results/metaannot.log   # the headline per dataset"
