#!/bin/bash
# run_us_prep.sh - stage 1: pull the system window-by-window (100 ps each),
# each window starting from the restart of the previous one. Serial, but short.
# Output: prep/umb_<centre>.rst = starting structure of each production window.
set -e
source ../settings.sh
source ./windows.sh
TEMPLATE=${TEMPLATE:-plumed_us.template}   # TASS reuses this script (see Section 6)
mkdir -p prep && cd prep

run_branch () {
  local prev=$1; shift
  for c in "$@"; do
    sed "s/XXXX/$c/" ../$TEMPLATE > plumed.dat
    $MDEXE -O -i ../us_prep.in -o umb_$c.out -p $PRMTOP -c $prev -r umb_$c.rst
    mv COLVAR COLVAR_prep_$c
    prev=umb_$c.rst
    echo "prepared window $c"
  done
}
run_branch $START_RST $BRANCH1
run_branch umb_-3.2.rst $BRANCH2
