# lib.sh - helpers sourced by every job script.
_LIBDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$_LIBDIR/../config.sh"
source "$_LIBDIR/windows.sh"

# hours available to md.py inside a job = SB_WALLTIME - STOP_MARGIN_MIN
walltime_hours () {
  local w=${1:-$SB_WALLTIME} d=0 h m s
  if [[ $w == *-* ]]; then d=${w%%-*}; w=${w#*-}; fi
  IFS=: read -r h m s <<< "$w"
  python3 -c "print(max(0.1, $d*24 + $h + $m/60 + ${s:-0}/3600 - $STOP_MARGIN_MIN/60))"
}

# make_plumed TEMPLATE OUTFILE [CENTRE]
#   fills @CENTRE@, @CHK_STEPS@ and @OPES_LOAD@ (LOAD line for the OPES library)
make_plumed () {
  local tpl=$1 out=$2 c=${3:-}
  local chk_steps=$(python3 -c "print(int(round($CHK_PS*1000/$DT_FS)))")
  local load=""
  if [ -n "$OPES_LIB" ]; then
    [ -f "$OPES_LIB" ] || { echo "OPES library $OPES_LIB missing: run common/build_opes.sh" >&2; return 1; }
    load="LOAD FILE=$OPES_LIB"
  fi
  sed -e "s#@CENTRE@#$c#g" -e "s#@CHK_STEPS@#$chk_steps#g" -e "s#@OPES_LOAD@#$load#g" "$tpl" > "$out"
}

# run_md OUTDIR PLUMED NS START  -> md.py with the job's wall-time budget; returns md.py's exit code
run_md () {
  python3 -u "$CODE/common/md.py" --prmtop "$PRMTOP" --outdir "$1" --plumed "$2" \
      --ns "$3" --start "$4" --max-hours "$(hours_left)"
}

# call once at the start of a job: fixes the deadline for all md.py runs in the job
start_clock () {
  export DEADLINE=$(python3 -c "import time; print(int(time.time() + 3600*$(walltime_hours)))")
}

# wall time left in this job (h)
hours_left () {
  if [ -n "${MAX_HOURS:-}" ]; then echo "$MAX_HOURS"; return; fi
  if [ -z "${DEADLINE:-}" ]; then walltime_hours; return; fi
  python3 -c "import time; print(max(0.0, ($DEADLINE - time.time())/3600))"
}

# resubmission guard: at most MAX_RESUBMIT continuation jobs per output folder
may_resubmit () {
  local f="$1/resubmits" n=0
  [ -f "$f" ] && n=$(cat "$f")
  if [ "$n" -ge "${MAX_RESUBMIT:-100}" ]; then
    echo "resubmission limit reached for $1"; return 1
  fi
  echo $((n + 1)) > "$f"
}
