#!/usr/bin/env bash
# CLAUDE.md rule 1: run doctor before any long run, and READ its output.
# CLAUDE.md rule 2: --install-plan, never --fix. The plan is a script to read BEFORE
# anything touches the filesystem. This script downloads nothing.
set -uo pipefail
cd "$(dirname "$0")"
MA=${MA:-$HOME/metaannot/metaannot.py}
PY=${PY:-python3}
rc=0
for d in */; do d=${d%/}
  [ -f "$d/config.yaml" ] || continue
  echo "================ $d ================"
  "$PY" "$MA" doctor --config "$d/config.yaml" --install-plan "$d/install.sh" 2>&1 | tee "$d/doctor.log"
  s=${PIPESTATUS[0]}; [ "$s" -ne 0 ] && rc=1
done
echo
echo "=================================================================="
[ "$rc" -eq 0 ] && echo "every dataset satisfied. Go to 02_run.sh." || cat <<'MSG'
Something is missing. Before installing anything:
  * The db: block of each config marks four paths INFERRED (pfam_hmm, dbcan_hmm,
    kofam_profiles/ko_list, interproscan_sh). If doctor says MISS for one of those,
    the database is probably present under another name — fix the PATH in the config
    first, and only download if it is genuinely absent.
  * Read <dataset>/install.sh before running it. It prints two totals, to download and
    on disk, and the per-item sizes are hand-maintained: check anything over ~10 GB
    against the provider.
  * MANUAL items (SignalP 6.0, InterProScan) are never attempted by doctor.
MSG
exit $rc
