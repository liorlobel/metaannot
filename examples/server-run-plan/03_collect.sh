#!/usr/bin/env bash
# One table across every finished dataset, and the KOfam control that motivated the run.
set -uo pipefail
cd "$(dirname "$0")"
printf "%-38s %10s %10s %9s %8s %8s\n" dataset 1_ko_path 2_ko_orph 3_ann_noKO 3d_duf 4_dark
for d in */; do d=${d%/}; s="$d/results/bin_summary.tsv"; [ -f "$s" ] || continue
  awk -F'\t' -v n="$d" '$1=="1_ko_pathway"{a=$3} $1=="2_ko_orphan"{b=$3} $1=="3_annotated_no_ko"{c=$3}
    $1=="3d_duf_only"{e=$3} $1=="4_dark"{f=$3}
    END{printf "%-38s %9s%% %9s%% %8s%% %7s%% %7s%%\n", n,a,b,c,e,f}' "$s"
done
echo
echo "== KOfam vs eggNOG: the control =="
grep -h "KOfamScan:" */results/metaannot.log 2>/dev/null || echo "  (no kofam lines — did the stage run?)"
echo
echo "== KEGG-invisible fraction =="
for d in */; do d=${d%/}; l="$d/results/metaannot.log"; [ -f "$l" ] || continue
  printf "  %-38s %s\n" "$d" "$(grep -oE '[0-9.]+% of quantified protein groups are invisible' "$l" | tail -1)"
done
echo
echo "== anything that must not be ignored =="
grep -l "parsed 0" */results/metaannot.log 2>/dev/null && echo "  ^ 'parsed 0 from a non-empty file' is NEVER noise: a tool's output was unreadable and every protein silently lost that evidence. Investigate before using these numbers."
grep -h "SEVER\|FATAL" */results/metaannot.log 2>/dev/null | sort -u | head
