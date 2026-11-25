#!/bin/bash
# job_check.sh - run on a GPU node BEFORE any production:
#   1. checks/check_env.py   (software, OPES, atoms, start state, PLUMED forces OpenCL vs CPU)
#   2. smoke + restart test of every method: 10 ps, forced stop, restart, 10 ps more
#      (checkpoints every 5 ps), then checks/check_outputs.py on the results
#   3. GPU speed of each method -> estimated days for the full production
# Results in $PROJ/00_checks and in the job log. Nothing here is used by production.
set -uo pipefail
source "${CODE:?}/common/lib.sh"
activate_env
echo "job ${SLURM_JOB_ID:-local} on $(hostname), GPU(s): ${CUDA_VISIBLE_DEVICES:-?}"
nvidia-smi --query-gpu=name,driver_version,memory.used --format=csv 2>/dev/null || true

python3 "$CODE/checks/check_env.py"; env_rc=$?

SMK=$PROJ/00_checks; rm -rf "$SMK"; mkdir -p "$SMK"
export CHK_PS=5 REPORT_PS=1 DCD_PS=0          # short intervals for the test only
for m in us metad tamd tass; do
  d=$SMK/$m; mkdir -p "$d"
  case $m in
    us)   make_plumed "$CODE/common/plumed/us.template"   "$d/plumed.dat" -1.4 ;;
    tass) make_plumed "$CODE/common/plumed/tass.template" "$d/plumed.dat" -1.4 ;;
    *)    make_plumed "$CODE/common/plumed/$m.dat"        "$d/plumed.dat" ;;
  esac || { echo "[FAIL] $m: PLUMED input"; continue; }
  # OPES/MetaD deposit every 500 steps = 1 ps, so 20 ps gives ~20 kernels/hills
  for part in 1 2; do
    extra=""; [ $part -eq 1 ] && extra="--stop-after-chunks 2"
    python3 -u "$CODE/common/md.py" --prmtop "$PRMTOP" --outdir "$d" --plumed "$d/plumed.dat" \
        --ns 0.02 --start "$START_XML" --max-hours 0.9 $extra > "$d/part$part.log" 2>&1
    rc=$?
    want=$([ $part -eq 1 ] && echo 3 || echo 0)
    [ $rc -eq $want ] && echo "[OK  ] $m part $part (exit $rc)" \
                      || { echo "[FAIL] $m part $part exit $rc (want $want): tail of log"; tail -25 "$d/part$part.log"; }
  done
done
python3 "$CODE/checks/check_outputs.py" "$SMK"; out_rc=$?
echo; echo "environment check exit $env_rc, output check exit $out_rc"
[ $env_rc -eq 0 ] && [ $out_rc -eq 0 ] && echo "ALL CHECKS PASSED" || echo "SOME CHECKS FAILED: do not start production"
