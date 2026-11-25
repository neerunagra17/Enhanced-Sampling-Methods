#!/bin/bash
# run_fes_checkpoints.sh - FES at all t_eff checkpoints for US, TASS, WT-MetaD and TAMD.
# Usage:  DATA=/path/to/production_output bash analysis/run_fes_checkpoints.sh
#   DATA must contain US/win_<c>/COLVAR, TASS/win_<c>/COLVAR, WT_MetaD/{HILLS,COLVAR}, TAMD/COLVAR.
# Output goes to results/ (run from anywhere).
set -e
cd "$(dirname "$0")/../results"
A=../analysis
D=${DATA:?set DATA to the folder with the production output}
export PYTHONPATH=$A
TEFF="10,20,40,60,80,100,150,200,300,500"
T_LIST=${TEFF//,/ }
source ../common/windows.sh
NWIN=$(echo $ALL | wc -w); TMIN=500
mkdir -p meta
: > meta/us.dat; : > meta/tass.dat
for c in $ALL; do
  echo "$D/US/win_$c/COLVAR $c $KAPPA_US"   >> meta/us.dat
  echo "$D/TASS/win_$c/COLVAR $c $KAPPA_US" >> meta/tass.dat
done
job() { m=$1; t=$2
  tw=$(python3 -c "print($t*1000/$NWIN - $TMIN)")
  python3 -c "import sys; sys.exit(0 if $tw > 0 else 1)" || { echo "$m $t ns: skip (too short for $NWIN windows)"; return; }
  if [ $m = us ]; then
    python3 $A/wham2d.py --meta meta/us.dat --T 300 --tmin $TMIN --tw $tw --out fes_us_${t}ns.dat > log_us_$t.txt
  else
    python3 $A/wham2d.py --meta meta/tass.dat --phi ext.phi_fict --psi ext.psi_fict --logw metad.rbias \
       --T 900 --tmin $TMIN --tw $tw --out fes_tass_${t}ns.dat > log_tass_$t.txt
  fi; echo "$m $t done"; }
export -f job; export A NWIN TMIN
( for t in $T_LIST; do echo "us $t"; echo "tass $t"; done | xargs -P 8 -L1 bash -c 'job $0 $1' ) &
python3 $A/sumhills_py.py --hills $D/WT_MetaD/HILLS --teff $TEFF --prefix fes_metad > log_metad.txt
python3 $A/reweight_fes.py --colvar $D/WT_MetaD/COLVAR --bias metad.rbias --T 300 --teff $TEFF --prefix fes_metadrw --tmin 1000 > log_metadrw.txt
python3 $A/reweight_fes.py --colvar $D/TAMD/COLVAR --phi ext.phi_fict --psi ext.psi_fict --T 1500 --teff $TEFF --prefix fes_tamd > log_tamd.txt
wait
echo ALL DONE
