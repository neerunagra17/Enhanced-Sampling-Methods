#!/bin/bash
# run_tamd.sh - one continuous TAMD run (100 ns, prod_tamd.in).
# Kept as a single run: check with the lab whether their EXTENDED_LAGRANGIAN
# can restart the auxiliary variables before splitting it into segments.
set -e
source ../settings.sh
cp plumed_tamd.dat plumed.dat
$MDEXE -O -i prod_tamd.in -o prod.out -p $PRMTOP -c $START_RST -r prod.rst -x prod.nc
