#!/bin/bash
# run_equil.sh - minimise, heat and equilibrate the solvated dipeptide.
# Run inside $PROJ/system after tleap (or skip entirely if the lab gives ala_npt.rst).
set -e
source ../settings.sh

$MDEXE -O -i min.in  -o min1.out -p ala_wat.prmtop -c ala_wat.inpcrd -r min1.rst -ref ala_wat.inpcrd
sed 's/ntr=1/ntr=0/' min.in > min2.in
$MDEXE -O -i min2.in -o min2.out -p ala_wat.prmtop -c min1.rst      -r min2.rst
$MDEXE -O -i heat.in -o heat.out -p ala_wat.prmtop -c min2.rst      -r heat.rst -x heat.nc -ref min2.rst
$MDEXE -O -i npt.in  -o npt.out  -p ala_wat.prmtop -c heat.rst      -r ala_npt.rst -x npt.nc

# quick checks: density should settle near 1.0 g/cm3, temperature near 300 K
grep -E "Density" npt.out | tail -3
grep -E "TEMP\(K\)" npt.out | tail -3
