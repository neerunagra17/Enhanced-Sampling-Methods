#!/bin/bash
# run_buckets.sh - one independent OPES run per bucket (folder bucket_<k>).
# Fill BUCKETS with "XI_L:XI_U" pairs and WALL_K from the group bucket set-up.
# Starting structures: one equilibrated structure inside each bucket
# (e.g. take prep/umb_<centre>.rst from the umbrella preparation).
set -e
source ../settings.sh
WALL_K=${WALL_K:-REPLACE_ME}
BUCKETS=("REPLACE_XI_L1:REPLACE_XI_U1" "REPLACE_XI_L2:REPLACE_XI_U2" "REPLACE_XI_L3:REPLACE_XI_U3")
STARTS=("../us/prep/umb_-2.0.rst" "../us/prep/umb_0.0.rst" "../us/prep/umb_2.0.rst")
for k in "${!BUCKETS[@]}"; do
  IFS=: read xl xu <<< "${BUCKETS[$k]}"
  d=bucket_$((k+1)); mkdir -p $d; cd $d
  sed -e "s/XI_L/$xl/" -e "s/XI_U/$xu/" -e "s/WALL_K/$WALL_K/g" ../plumed_bucket.template > plumed.dat
  $MDEXE -O -i ../bucket_prod.in -o prod.out -p $PRMTOP -c ../${STARTS[$k]} -r prod.rst -x prod.nc
  cd ..
done
