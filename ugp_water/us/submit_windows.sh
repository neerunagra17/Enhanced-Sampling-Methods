#!/bin/bash
# submit_windows.sh - submit every production window as its own job (parallel).
# Run inside us/ (or tass/ with TEMPLATE/MDIN set, see Section 6).
source ./windows.sh
TEMPLATE=${TEMPLATE:-plumed_us.template}
MDIN=${MDIN:-us_prod.in}
for c in $ALL; do
  d=win_$c; mkdir -p $d
  sed "s/XXXX/$c/" $TEMPLATE > $d/plumed.dat
  CMD="cd $PWD/$d && \$MDEXE -O -i ../$MDIN -o prod.out -p \$PRMTOP -c ../prep/umb_$c.rst -r prod.rst -x prod.nc"
  qsub -N w$c -v CMD="$CMD" ../job.pbs
done
