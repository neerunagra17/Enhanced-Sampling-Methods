#!/bin/bash
# job_single.sh metad|tamd - one 500 ns trajectory from the equilibrated NPT state.
# Resubmits itself until finished; every continuation restarts from the last checkpoint.
set -uo pipefail
source "${CODE:?}/common/lib.sh"
activate_env
start_clock
M=$1
case $M in
  metad) DIR=$METAD_DIR; TPL=$CODE/common/plumed/metad.dat ;;
  tamd)  DIR=$TAMD_DIR;  TPL=$CODE/common/plumed/tamd.dat ;;
  *) echo "metad|tamd"; exit 1 ;;
esac
mkdir -p "$DIR"
echo "job ${SLURM_JOB_ID:-local} on $(hostname), GPU(s): ${CUDA_VISIBLE_DEVICES:-?}; $M -> $DIR"
# the PLUMED input is written once; later jobs must use the identical file
[ -f "$DIR/plumed.dat" ] || make_plumed "$TPL" "$DIR/plumed.dat" || exit 1
run_md "$DIR" "$DIR/plumed.dat" "$SINGLE_NS" "$START_XML"
rc=$?
if [ $rc -eq 3 ]; then
  may_resubmit "$DIR" && bash "$CODE/submit.sh" $M
elif [ $rc -ne 0 ]; then
  echo "$M failed (exit $rc)"; exit $rc
fi
