#!/bin/bash
# job_prep.sh us|tass - chained 100 ps preparation of the 32 windows:
# window c is restrained at phi = c (physical phi, for US and TASS alike) and starts from the
# final state of the previous window; branch 1 starts from the equilibrated NPT state,
# branch 2 from the -3.2 window. Finished windows (DONE) are skipped, so the job can resume.
# Preparation time (32 x 100 ps = 3.2 ns) is NOT part of t_eff; it is reported separately.
set -uo pipefail
source "${CODE:?}/common/lib.sh"
activate_env
start_clock
M=$1
case $M in us) DIR=$US_DIR ;; tass) DIR=$TASS_DIR ;; *) echo "us|tass"; exit 1 ;; esac
P=$DIR/prep; mkdir -p "$P"
echo "job ${SLURM_JOB_ID:-local} on $(hostname), GPU(s): ${CUDA_VISIBLE_DEVICES:-?}; $M preparation in $P"
NS=$(python3 -c "print($PREP_PS/1000)")

prep_branch () {
  local prev=$1; shift
  for c in "$@"; do
    local d=$P/win_$c
    if [ ! -f "$d/DONE" ]; then
      mkdir -p "$d"
      make_plumed "$CODE/common/plumed/us.template" "$d/plumed.dat" "$c" || exit 1
      run_md "$d" "$d/plumed.dat" "$NS" "$prev"
      rc=$?
      if [ $rc -eq 3 ]; then
        may_resubmit "$P" && bash "$CODE/submit.sh" ${M}_prep
        exit 0
      elif [ $rc -ne 0 ]; then
        echo "window $c failed (exit $rc)"; exit $rc
      fi
    fi
    prev=$d/state.xml
    echo "prepared $c: final phi = $(tail -1 $d/COLVAR | awk '{print $2}')"
  done
}
prep_branch "$START_XML" $BRANCH1
prep_branch "$P/win_$BRANCH2_FROM/state.xml" $BRANCH2
echo "$M preparation complete: next  bash submit.sh ${M}_prod"
