#!/usr/bin/env bash
# Build the identified-protein FASTA for each dataset.
# Annotate what was identified, not the whole catalogue: it is far faster and the bin
# percentages then describe the actual result rather than the database.
# The UHGP-90 FASTA is 9.6 GB and iMGMC 2.4 GB, so this reads ~46 GB in total. Minutes, not hours.
set -euo pipefail
cd "$(dirname "$0")"
MA=${MA:-$HOME/metaannot/metaannot.py}
PY=${PY:-python3}
for d in */; do d=${d%/}
  [ -f "$d/config.yaml" ] || continue
  fa=$(cat "$d/FASTA_PATH")
  q=$(grep -m1 '^quant_table:' "$d/config.yaml" | cut -d'"' -f2)
  if [ -s "$d/input/proteins.faa" ]; then echo "== $d: already subset, skipping"; continue; fi
  echo "== $d"
  "$PY" "$MA" subset --db "$fa" --quant "$q" --format fragpipe_peptide \
       --out "$d/input/proteins.faa" 2>&1 | tee "$d/subset.log" | tail -3
done
echo
echo "sequence counts:"
for d in */; do d=${d%/}; [ -s "$d/input/proteins.faa" ] && printf "  %-38s %s\n" "$d" "$(grep -c '^>' "$d/input/proteins.faa")"; done
