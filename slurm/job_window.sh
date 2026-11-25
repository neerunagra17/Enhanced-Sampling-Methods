#!/bin/bash
# job_window.sh us|tass - one production window per array task (index 0..31 -> centre).
# Starts from prep/win_<c>/state.xml (positions, velocities, box of the end of the 100 ps
# preparation); the production clock starts at 0. Length WIN_NS (15.625 ns).
set -uo pipefail
source "${CODE:?}/common/lib.sh"
activate_env
start_clock
M=$1
case $M in
  us)   DIR=$US_DIR;   TPL=$CODE/common/plumed/us.template ;;
  tass) DIR=$TASS_DIR; TPL=$CODE/common/plumed/tass.template ;;
  *) echo "us|tass"; exit 1 ;;
esac
IDX=${SLURM_ARRAY_TASK_ID:?array task id missing}
c=$(centre_of $IDX)
OUT=$DIR/win_$c; mkdir -p "$OUT"
START=$DIR/prep/win_$c/state.xml
[ -f "$DIR/prep/win_$c/DONE" ] || { echo "preparation of window $c not finished"; exit 1; }
echo "job ${SLURM_JOB_ID:-local} task $IDX on $(hostname), GPU(s): ${CUDA_VISIBLE_DEVICES:-?}; $M window $c -> $OUT"
[ -f "$OUT/plumed.dat" ] || make_plumed "$TPL" "$OUT/plumed.dat" "$c" || exit 1
run_md "$OUT" "$OUT/plumed.dat" "$WIN_NS" "$START"
rc=$?
if [ $rc -eq 3 ]; then
  may_resubmit "$OUT" && bash "$CODE/submit.sh" ${M}_prod $IDX
elif [ $rc -ne 0 ]; then
  echo "window $c failed (exit $rc)"; exit $rc
fi
