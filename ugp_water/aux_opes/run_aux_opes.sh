#!/bin/bash
# run_aux_opes.sh - OPTIONAL auxiliary-CV + OPES run (single 100 ns run).
set -e
source ../settings.sh
cp plumed_aux_opes.dat plumed.dat
$MDEXE -O -i prod_aux_opes.in -o prod.out -p $PRMTOP -c $START_RST -r prod.rst -x prod.nc
