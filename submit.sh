#!/bin/bash
# submit.sh - submit every job of the project with the SLURM settings in config.sh.
#
#   bash submit.sh check                 GPU environment check + 10 ps smoke/restart test of every method
#   bash submit.sh us_prep | tass_prep   chained 100 ps preparation of the 32 windows (one GPU job)
#   bash submit.sh us_prod | tass_prod [ARRAY]   production windows (job array, default 0-31)
#   bash submit.sh metad | tamd          500 ns single trajectory
#   bash submit.sh status                progress of every run (no job)
#
# Jobs that reach their wall-time limit stop at a checkpoint and resubmit themselves.
set -euo pipefail
export CODE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$CODE/common/lib.sh"
mkdir -p "$PROJ/logs"
target=${1:?usage: bash submit.sh TARGET [ARRAY]}

GPU="--partition=$SB_PARTITION --gres=$SB_GRES --mem=$SB_MEM --cpus-per-task=$SB_CPUS --ntasks=1 --time=$SB_WALLTIME --export=ALL,CODE=$CODE"
CPU="--partition=$SB_CPU_PARTITION --mem=$SB_CPU_MEM --cpus-per-task=2 --ntasks=1 --time=$SB_CPU_WALLTIME --export=ALL,CODE=$CODE"
[ -n "$SB_ACCOUNT" ] && { GPU="$GPU --account=$SB_ACCOUNT"; CPU="$CPU --account=$SB_ACCOUNT"; }
L=$PROJ/logs

prep_done () {   # all 32 preparation windows of $1 (us|tass) finished?
  local dir=$2 n=0
  for c in $ALL; do [ -f "$dir/prep/win_$c/DONE" ] && n=$((n+1)); done
  [ "$n" -eq "$NWIN" ] || { echo "$1 preparation incomplete ($n/$NWIN windows): run 'bash submit.sh $1_prep' first"; return 1; }
}

case $target in
  check)
    sbatch $GPU --time=01:00:00 -J ugp_check -o $L/check_%j.out "$CODE/slurm/job_check.sh" ;;
  us_prep)
    sbatch $GPU -J us_prep -o $L/us_prep_%j.out "$CODE/slurm/job_prep.sh" us ;;
  tass_prep)
    sbatch $GPU -J tass_prep -o $L/tass_prep_%j.out "$CODE/slurm/job_prep.sh" tass ;;
  us_prod)
    prep_done us "$US_DIR"
    sbatch $GPU --array=${2:-0-$((NWIN-1))} -J us_win -o $L/us_win_%a_%j.out "$CODE/slurm/job_window.sh" us ;;
  tass_prod)
    prep_done tass "$TASS_DIR"
    sbatch $GPU --array=${2:-0-$((NWIN-1))} -J tass_win -o $L/tass_win_%a_%j.out "$CODE/slurm/job_window.sh" tass ;;
  metad|tamd)
    sbatch $GPU -J $target -o $L/${target}_%j.out "$CODE/slurm/job_single.sh" $target ;;
  status)
    python3 "$CODE/common/status.py" ;;
  *) echo "unknown target $target"; exit 1 ;;
esac
