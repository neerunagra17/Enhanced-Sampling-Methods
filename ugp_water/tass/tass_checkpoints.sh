#!/bin/bash
# tass_checkpoints.sh - TASS FES at the t_eff checkpoints (cross-check route b).
set -e
source ../settings.sh
source ./windows.sh
TMIN=${TMIN:-500}; TS=900
NWIN=$(echo $ALL | wc -w)
: > meta.dat
for c in $ALL; do echo "win_$c/COLVAR  $c  $KAPPA" >> meta.dat; done
for t in ${TEFF//,/ }; do
  tw=$(python3 -c "print($t*1000/$NWIN - $TMIN)")
  python3 -c "import sys; sys.exit(0 if $tw > 0 else 1)" || { echo "skip $t ns"; continue; }
  python3 $ANALYSIS/wham2d.py --meta meta.dat --phi ext.phi_fict --psi ext.psi_fict \
      --logw metad.rbias --T $TS --nbins $NBINS --tmin $TMIN --tw $tw \
      --out ../results/fes_tass_${t}ns.dat
done
