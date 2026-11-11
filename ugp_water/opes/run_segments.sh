#!/bin/bash
# run_segments.sh - run a single-trajectory method as NSEG segments of 10 ns.
# Usage (inside the method folder):  bash run_segments.sh NSEG FIRST_PLUMED RESTART_PLUMED
#   e.g. bash run_segments.sh 20 plumed_metad.dat plumed_metad_restart.dat   (200 ns)
# Segment 1 starts from the equilibrated structure, later segments from the
# previous restart, with the RESTART version of plumed.dat so the bias continues.
# Already finished segments (seg_<k>.rst exists) are skipped, so the script can
# simply be resubmitted after a cluster time-out.
set -e
source ../settings.sh
NSEG=$1; FIRST=$2; RESTARTFILE=$3
MDIN=${MDIN:-prod_segment.in}
prev=$START_RST
for k in $(seq 1 $NSEG); do
  if [ -f seg_$k.rst ]; then prev=seg_$k.rst; continue; fi
  if [ $k -eq 1 ]; then cp $FIRST plumed.dat; else cp $RESTARTFILE plumed.dat; fi
  $MDEXE -O -i $MDIN -o seg_$k.out -p $PRMTOP -c $prev -r seg_$k.rst -x seg_$k.nc
  prev=seg_$k.rst
  echo "segment $k done"
done
