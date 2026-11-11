#!/bin/bash
# run_us_prod.sh - stage 2: production of all windows, INDEPENDENT of each other.
# Each window runs in its own folder win_<centre>/ (PLUMED writes COLVAR in the
# working directory, so windows must not share a folder).
# On a cluster, submit one job per window instead of the loop (see Tutorial).
set -e
source ../settings.sh
source ./windows.sh
TEMPLATE=${TEMPLATE:-plumed_us.template}
MDIN=${MDIN:-us_prod.in}

for c in $ALL; do
  d=win_$c; mkdir -p $d; cd $d
  sed "s/XXXX/$c/" ../$TEMPLATE > plumed.dat
  $MDEXE -O -i ../$MDIN -o prod.out -p $PRMTOP -c ../prep/umb_$c.rst \
         -r prod.rst -x prod.nc
  cd ..
  echo "finished window $c"
done
