#!/bin/bash
# run_l2_all.sh - L2 error vs t_eff for every method, against
#   (a) the independent reference  ref/fes_ref.dat
#   (b) each method's own longest run (internal convergence)
# and the two summary plots for the report. Run inside $PROJ/results.
set -e
source ../settings.sh
METHODS="us metad opes tamd tass opes_bucket"
last=$(echo ${TEFF//,/ } | awk '{print $NF}')
ref_csv=""; int_csv=""
for m in $METHODS; do
  files=$(ls fes_${m}_*ns.dat 2>/dev/null | grep -E "^fes_${m}_[0-9.]+ns\.dat$" || true)
  [ -z "$files" ] && { echo "no files for $m, skipping"; continue; }
  M=$(echo $m | tr a-z A-Z | sed 's/OPES_BUCKET/OPES+bucket/')
  python3 $ANALYSIS/l2_error.py --ref ../ref/fes_ref.dat --sim $files \
      --fcut $FCUT --nbins $NBINS --method "$M" --out l2_${m}_ref.csv
  ref_csv="$ref_csv l2_${m}_ref.csv"
  if [ -f fes_${m}_${last}ns.dat ]; then
    python3 $ANALYSIS/l2_error.py --ref fes_${m}_${last}ns.dat --sim $files \
        --fcut $FCUT --nbins $NBINS --method "$M" --out l2_${m}_int.csv
    int_csv="$int_csv l2_${m}_int.csv"
  fi
done
python3 $ANALYSIS/plot_l2.py --csv $ref_csv --title "Error w.r.t. reference" --out l2_vs_teff_ref.pdf
[ -n "$int_csv" ] && python3 $ANALYSIS/plot_l2.py --csv $int_csv \
    --title "Internal convergence" --out l2_vs_teff_internal.pdf
