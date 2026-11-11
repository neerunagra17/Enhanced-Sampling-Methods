# settings.sh - shared settings, sourced by every run/analysis script.
# Edit ONCE here so all methods use identical choices.

export PROJ=$HOME/ugp_water                 # project root (see project_layout.txt)
export ANALYSIS=$PROJ/analysis
export PRMTOP=$PROJ/system/ala_wat.prmtop   # topology (solvated)
export START_RST=$PROJ/system/ala_npt.rst   # equilibrated starting structure

# MD engine: sander.MPI works with PLUMED in AMBER >= 20 (plumed=1 in &cntrl).
# Use pmemd only if your lab's build is known to support PLUMED.
export NPROC=8
export MDEXE="mpirun -np $NPROC $AMBERHOME/bin/sander.MPI"

# common analysis choices (energies in kJ/mol, angles in rad)
export TEMP=300                     # physical temperature, K
export NBINS=72                     # FES grid 72 x 72 (5 degrees)
export FCUT=40                      # kJ/mol, cut-off in the L2 error (fixed with supervisor)
export TEFF="10,20,40,60,80,100,150,200"   # t_eff checkpoints, ns (fixed with supervisor)
