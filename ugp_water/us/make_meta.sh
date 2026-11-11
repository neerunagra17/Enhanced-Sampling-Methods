#!/bin/bash
# make_meta.sh - write meta.dat for wham2d.py / plot_overlap.py:
#   path_to_COLVAR   centre(rad)   kappa(kJ/mol/rad^2)
source ./windows.sh
: > meta.dat
for c in $ALL; do echo "win_$c/COLVAR  $c  $KAPPA" >> meta.dat; done
echo "wrote meta.dat ($(wc -l < meta.dat) windows)"
