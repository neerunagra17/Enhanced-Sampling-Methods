#!/bin/bash
# us_checkpoints.sh - US free-energy surfaces at the common t_eff checkpoints.
#   t_eff = (number of windows) x (time simulated per window)
# For a checkpoint t_eff, each window contributes its first t_eff/NWIN of
# production; the first TMIN ps of that are discarded as equilibration
# (they still count in t_eff, so all methods are charged for all MD time).
set -e
source ../settings.sh
source ./windows.sh
TMIN=${TMIN:-500}                       # ps discarded per window
NWIN=$(echo $ALL | wc -w)
mkdir -p ../results
for t in ${TEFF//,/ }; do
  tw=$(python3 -c "print($t*1000/$NWIN - $TMIN)")     # ps of data used per window
  if python3 -c "import sys; sys.exit(0 if $tw > 0 else 1)"; then
    python3 $ANALYSIS/wham2d.py --meta meta.dat --T $TEMP --nbins $NBINS \
        --tmin $TMIN --tw $tw --out ../results/fes_us_${t}ns.dat
  else
    echo "t_eff = $t ns too short for $NWIN windows (needs > $TMIN ps each)"
  fi
done
