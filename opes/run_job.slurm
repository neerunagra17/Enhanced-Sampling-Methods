#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
# GPU launcher for OPES runs with openmm_run.py (run from the job folder):
#   sbatch --job-name=NAME --export=ALL,NS=500 run_job.slurm
# Needs in the folder: topol.top, start.gro, plumed_first.dat, plumed_restart.dat, openmm_run.py, test_ctx.py
# Tries CUDA -> CUDA with NVRTC 12.2 (old-driver fix) -> OpenCL, uses the first that creates a context.
# Set MM_ROOT (micromamba installation with an environment "omm": OpenMM, openmm-plumed, PLUMED 2.10)
# and, if needed, NVRTC_DIR (NVRTC 12.2 libraries) before submitting.
cd "$SLURM_SUBMIT_DIR"
: "${MM_ROOT:?set MM_ROOT to the micromamba installation}"
export MAMBA_ROOT_PREFIX=$MM_ROOT/root
eval "$($MM_ROOT/bin/micromamba shell hook -s bash)"
micromamba activate omm
unset PLUMED_ROOT
export PLUMED_KERNEL=$CONDA_PREFIX/lib/libplumedKernel.so
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
echo "host $(hostname)  job $SLURM_JOB_ID  $(date)"; nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
NV=${NVRTC_DIR:-$MM_ROOT/nvrtc122/lib}
PLAT=""
for mode in A B C; do
  case $mode in
    A) P=CUDA; unset LD_PRELOAD ;;
    B) P=CUDA; export LD_PRELOAD=$NV/libnvrtc.so.12; export LD_LIBRARY_PATH=$NV:$LD_LIBRARY_PATH ;;
    C) P=OpenCL; unset LD_PRELOAD ;;
  esac
  if python test_ctx.py $P 2>&1 | tail -3 | tee /dev/stderr | grep -q "context OK"; then PLAT=$P; echo "using mode $mode ($P)"; break; fi
done
[ -z "$PLAT" ] && { echo "no working GPU platform"; exit 1; }
python openmm_run.py --format gromacs --top topol.top --coords start.gro --ns ${NS:-500} --chunk-ns ${CHUNK:-2} \
       --max-hours ${MAXH:-44} --platform $PLAT --dt-fs ${DT:-2} ${HMR:+--hmr $HMR}
rc=$?
if [ $rc -eq 3 ]; then
  n=$(cat resubmits 2>/dev/null || echo 0)
  if [ "$n" -lt 6 ]; then echo $((n+1)) > resubmits; sbatch --export=ALL run_job.slurm; fi
fi
exit $rc
