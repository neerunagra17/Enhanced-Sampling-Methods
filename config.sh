# config.sh - ONE place for every setting. Sourced by all shell scripts.
# Edit only the block marked "SITE".

# ----------------------------------------------------------------- SITE
export PROJ=${PROJ:?set PROJ to the project root, e.g. export PROJ=/scratch/<user>/alanine_water}
export CODE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)   # this folder (found automatically)
export CONDA_ENV=ugp_gpu
export PLATFORM=${PLATFORM:-OpenCL}                     # OpenCL (CUDA has driver problems); CPU for tests
export PRECISION=mixed                                  # OpenCL/CUDA precision
# activate the environment inside batch jobs (adjust if conda is a module on the cluster)
activate_env () {
  [ -n "${UGP_NO_ACTIVATE:-}" ] && return 0   # testing outside SLURM with the env already active
  # module load miniconda3       # <- uncomment / change if conda is provided as a module
  if [ -z "${CONDA_EXE:-}" ] && [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -n "${CONDA_EXE:-}" ]; then
    source "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
  fi
  conda activate $CONDA_ENV
}
# SLURM: values given for GPU jobs on this cluster. Check the partition time limit with
#   sinfo -p gpu -o "%P %l"     and set WALLTIME below it.
export SB_PARTITION=gpu
export SB_GRES=gpu:1
export SB_MEM=4G
export SB_CPUS=2              # PLUMED runs on the CPU every step; 2 cores is enough
export SB_WALLTIME=${SB_WALLTIME:-24:00:00}
export SB_ACCOUNT=""          # set if the cluster requires  #SBATCH -A <account>
export STOP_MARGIN_MIN=20     # stop cleanly this many minutes before the wall-time limit

# ----------------------------------------------------------------- SYSTEM (already prepared)
export SYS=$PROJ/01_system
export PRMTOP=$SYS/ala_water.prmtop
export START_XML=$SYS/gpu_npt_final.xml        # equilibrated NPT state (positions, velocities, box)

# ----------------------------------------------------------------- MD settings (must match the equilibration)
export TEMP=300               # K
export DT_FS=2                # time step, fs
export FRICTION=1.0           # Langevin friction, 1/ps
export CUTOFF_NM=1.0          # PME real-space cut-off, nm (10 A)
export EWALD_TOL=0.0005       # OpenMM default
export CHK_PS=${CHK_PS:-100}             # checkpoint (state + PLUMED sync point) every 100 ps
export REPORT_PS=${REPORT_PS:-10}           # log energies / speed every 10 ps
export DCD_PS=${DCD_PS:-10}              # trajectory frame every 10 ps (0 = no trajectory)

# ----------------------------------------------------------------- output folders (separate per method)
export US_DIR=$PROJ/02_us
export METAD_DIR=$PROJ/03_metad
export TAMD_DIR=$PROJ/05_tamd
export TASS_DIR=$PROJ/06_tass
export RES_DIR=$PROJ/07_results
export ANALYSIS=$CODE/analysis

# ----------------------------------------------------------------- production lengths (t_eff = 500 ns per method)
export NWIN=32
export PREP_PS=${PREP_PS:-100}            # chained preparation per window (not counted in t_eff, reported separately)
export WIN_NS=${WIN_NS:-15.625}          # production per window: 32 x 15.625 = 500 ns
export SINGLE_NS=${SINGLE_NS:-500}          # WT-MetaD, TAMD: one 500 ns trajectory each

# ----------------------------------------------------------------- analysis (energies kJ/mol, angles rad)
export NBINS=${NBINS:-72}               # common 72 x 72 grid (5 degrees)
export FCUT=${FCUT:-61.0}              # kJ/mol (= 14.58 kcal/mol), L2 cut-off on F_ref
export TEFF=${TEFF:-"10,20,40,60,80,100,150,200,300,500"}     # ns; checkpoints for every method
export TMIN_WIN_PS=${TMIN_WIN_PS:-500}        # ps discarded at the start of each US/TASS window (still counted in t_eff)

# ----------------------------------------------------------------- PLUMED
# conda-forge plumed 2.9.2 is built WITHOUT the opes module; common/build_opes.sh builds it as
# a library and OPES inputs then start with this LOAD line. Set to "" if your PLUMED has opes.
export OPES_LIB=$CODE/common/plumed_lib/OPESmetad.so

# ----------------------------------------------------------------- CPU jobs (analysis)
export SB_CPU_PARTITION=${SB_CPU_PARTITION:-standard}   # check with `sinfo`; analysis needs no GPU
export SB_CPU_MEM=8G
export SB_CPU_WALLTIME=12:00:00

# ----------------------------------------------------------------- reference FES for the L2 error
# 3 columns "phi psi F" (rad, rad, kJ/mol), same force field / water / 300 K, independent and
# much longer than 500 ns. NOT produced by this code: see README ("Reference surface").
export REF_FES=$PROJ/08_reference/fes_ref.dat
