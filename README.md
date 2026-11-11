# Enhanced Sampling of Alanine Dipeptide in Water: A Performance Comparison

**CHM 391 Undergraduate Project, Department of Chemistry, IIT Kanpur**
Kruti Marthak · Neeru Nagra  |  Supervisor: Prof. Nisanth N. Nair

This repository contains the input files, run scripts and analysis code used to compute the conformational free-energy surface (FES) $F(\phi,\psi)$ of alanine dipeptide in explicit water with six enhanced-sampling methods. It also contains the scripts that compare how quickly each method converges.

## Methods compared

| Method | Bias / mechanism | FES estimator | Folder |
|---|---|---|---|
| Umbrella sampling (US) | 32 harmonic windows along $\phi$ | 2D WHAM | `us/` |
| Well-tempered metadynamics (WT-MetaD) | Gaussian bias on $(\phi,\psi)$ | `sum_hills` / reweighting | `metad/` |
| OPES | On-the-fly probability-based bias on $(\phi,\psi)$ | Reweighting | `opes/` |
| TAMD / d-AFED | Auxiliary CVs at high temperature $T_s$ | Histogram at $T_s$ | `tamd/` |
| TASS | US on aux. $\phi$ + WT-MetaD on aux. $\psi$, at $T_s$ | Reweighting + WHAM | `tass/` |
| OPES + buckets | Flat-bottomed walls along $\phi$, OPES inside | Mean-force reconstruction | `opes_bucket/` |
| *Aux. CV + OPES (optional)* | OPES bias on TAMD auxiliary variables | Reweighting | `aux_opes/` |

**Performance metric:** the $L^2$ (RMS) error of each FES against a reference surface, plotted against the total simulation time $t_{\text{eff}} = n_{\text{windows}} \times t_{\text{window}}$:

$$
L^2 = \sqrt{\frac{1}{N_g}\sum_{i=1}^{N_g}\bigl[F_{\text{ref}}(i)-F_{\text{sim}}(i)\bigr]^2},\qquad F_{\text{ref}}(i) < F_{\text{cut}}
$$

All surfaces are put on a common 72 × 72 grid (5° bins) and shifted so that their minimum is zero.

## Simulation set-up

| | |
|---|---|
| System | Ace–Ala–Nme, AMBER ff14SB, TIP3P water, truncated octahedron |
| Engine | AMBER `sander.MPI` + PLUMED 2 |
| Ensemble | NVT, 300 K, Langevin thermostat (γ = 1 ps⁻¹) |
| Time step | 2 fs, SHAKE on bonds to hydrogen |
| Electrostatics | PME, 10 Å cut-off |
| CVs | $\phi$ (atoms 5-7-9-15), $\psi$ (atoms 7-9-15-17) |
| Units | kJ/mol, rad, ps throughout, as written by PLUMED |

## Repository layout

```
.
├── settings.sh        shared settings, sourced by every script (paths, NPROC, grid, F_cut, t_eff)
├── check_env.sh       checks AMBER, PLUMED (METAD, OPES_METAD, EXTENDED_LAGRANGIAN) and Python
├── job.pbs            generic PBS job wrapper
├── system/            tleap, minimisation, heating, NPT equilibration, atom-index check
├── us/                umbrella sampling: sequential preparation + parallel production windows
├── metad/             well-tempered metadynamics (restartable 10 ns segments)
├── opes/              OPES_METAD (restartable segments)
├── tamd/              TAMD / d-AFED with EXTENDED_LAGRANGIAN
├── tass/              TASS windows (aux. φ umbrellas + aux. ψ WT-MetaD at T_s = 900 K)
├── opes_bucket/       OPES with bucket sampling (template)
├── aux_opes/          optional: OPES on auxiliary variables
├── analysis/          Python analysis library
├── ref/               reference FES goes here (ref/fes_ref.dat)
└── results/           FES files, L2 tables and figures (run_l2_all.sh)
```

### Analysis scripts (`analysis/`)

| Script | Purpose |
|---|---|
| `fesutils.py` | Shared helpers: COLVAR reading, periodic grids, histograms, FES I/O, interpolation |
| `wham2d.py` | 2D WHAM for umbrella windows on a φ bias (also used for TASS at $T_s$) |
| `reweight_fes.py` | FES by reweighting (OPES, MetaD rbias, TAMD at $T_s$) at each $t_{\text{eff}}$ |
| `metad_checkpoints.py` | `plumed sum_hills` at each $t_{\text{eff}}$ checkpoint |
| `l2_error.py` | $L^2$ error of a set of FES files against a reference |
| `plot_l2.py` | $L^2$ vs $t_{\text{eff}}$ for all methods |
| `plot_fes.py`, `plot_overlap.py`, `plot_traces.py` | FES maps, window overlap, CV time series |

## Requirements

- AMBER (AmberTools ≥ 20) with `sander.MPI` patched with **PLUMED 2.8+**, including the `opes` module
- Python 3.8+ with `numpy`, `scipy` and `matplotlib`
- A PBS cluster (to use `job.pbs`); scripts also run interactively

```bash
pip install numpy scipy matplotlib
```

## Quick start

```bash
git clone https://github.com/neerunagra17/<repo-name>.git $HOME/ugp_water
cd $HOME/ugp_water
nano settings.sh          # set PROJ, NPROC, MDEXE
bash check_env.sh
```

**1. System:** use existing `system/ala_wat.prmtop` and `system/ala_npt.rst` files, or build them:
```bash
cd $PROJ/system
tleap -f tleap_water.in
qsub -v CMD="bash run_equil.sh" ../job.pbs
cpptraj -p ala_wat.prmtop -i check_atoms.cpptraj     # confirm φ/ψ atom indices
```

**2. Run each method:**
```bash
# Umbrella sampling
cd $PROJ/us && qsub -v CMD="bash run_us_prep.sh" ../job.pbs
bash submit_windows.sh                                   # after prep/ is complete

# WT-MetaD and OPES (20 × 10 ns segments, restartable)
cd $PROJ/metad && qsub -v CMD="bash run_segments.sh 20 plumed_metad.dat plumed_metad_restart.dat" ../job.pbs
cd $PROJ/opes  && qsub -v CMD="bash run_segments.sh 20 plumed_opes.dat plumed_opes_restart.dat" ../job.pbs

# TAMD
cd $PROJ/tamd && qsub -v CMD="bash run_tamd.sh" ../job.pbs

# TASS
cd $PROJ/tass && qsub -v CMD="bash ../us/run_us_prep.sh" ../job.pbs
TEMPLATE=plumed_tass.template MDIN=tass_prod.in bash ../us/submit_windows.sh

# OPES + buckets (set BUCKETS, WALL_K, STARTS first)
cd $PROJ/opes_bucket && qsub -v CMD="bash run_buckets.sh" ../job.pbs
```

**3. Free-energy surfaces at each $t_{\text{eff}}$:**
```bash
cd $PROJ/us    && bash make_meta.sh && bash us_checkpoints.sh
cd $PROJ/metad && python3 $ANALYSIS/metad_checkpoints.py --hills HILLS --teff $TEFF --prefix ../results/fes_metad
cd $PROJ/opes  && python3 $ANALYSIS/reweight_fes.py --colvar COLVAR --bias opes.bias --T 300 --teff $TEFF --prefix ../results/fes_opes
cd $PROJ/tamd  && python3 $ANALYSIS/reweight_fes.py --colvar COLVAR --phi ext.phi_fict --psi ext.psi_fict --T 1500 --teff $TEFF --prefix ../results/fes_tamd
cd $PROJ/tass  && bash tass_checkpoints.sh
```

**4. Performance comparison:** put the reference FES in `ref/fes_ref.dat` (columns `phi psi F`, kJ/mol), then:
```bash
cd $PROJ/results
bash run_l2_all.sh
```
This writes `l2_<method>_ref.csv` and `l2_<method>_int.csv`, and the two summary figures `l2_vs_teff_ref.pdf` (error against the reference) and `l2_vs_teff_internal.pdf` (internal convergence).

## Checks before trusting a surface

- **US / TASS:** neighbouring window histograms overlap over the full φ range (`plot_overlap.py`).
- **MetaD / OPES:** φ crosses between the φ < 0 and φ > 0 basins many times (`plot_traces.py`).
- **TAMD:** the physical CV follows the auxiliary variable closely (adiabatic separation).
- **Units:** PLUMED biases are in kJ/mol. Reweighting must use $k_BT$ in kJ/mol (2.494 at 300 K).
- **OPES:** `BARRIER` must exceed the φ barrier (≈ 40 kJ/mol). A value that is too low gives a stable but wrong FES.

## Status

- [x] Input files and analysis pipeline; analysis scripts validated on synthetic data with a known FES
- [x] Preliminary vacuum runs (US, TAMD, OPES)
- [ ] Production runs in water
- [ ] $L^2$ vs $t_{\text{eff}}$ comparison and final report

## Key references

1. G. M. Torrie and J. P. Valleau, *J. Comput. Phys.* **23**, 187 (1977). Umbrella sampling
2. S. Kumar *et al.*, *J. Comput. Chem.* **13**, 1011 (1992). WHAM
3. A. Barducci, G. Bussi and M. Parrinello, *Phys. Rev. Lett.* **100**, 020603 (2008). Well-tempered metadynamics
4. M. Invernizzi and M. Parrinello, *J. Phys. Chem. Lett.* **11**, 2731 (2020). OPES
5. L. Maragliano and E. Vanden-Eijnden, *Chem. Phys. Lett.* **426**, 168 (2006). TAMD
6. J. B. Abrams and M. E. Tuckerman, *J. Phys. Chem. B* **112**, 15742 (2008). d-AFED
7. S. Awasthi and N. N. Nair, *J. Chem. Phys.* **146**, 094108 (2017). TASS
8. A. Pal *et al.*, *J. Comput. Chem.* **42**, 1996 (2021). Mean-force TASS
9. R. Javed, A. B. Kapakayala and N. N. Nair, *J. Chem. Theory Comput.* **20** (2024), [doi:10.1021/acs.jctc.4c00776](https://pubs.acs.org/doi/10.1021/acs.jctc.4c00776). Bucket sampling
10. G. A. Tribello *et al.*, *Comput. Phys. Commun.* **185**, 604 (2014). PLUMED 2

The tutorials of the [NNN group, IIT Kanpur](https://sites.google.com/view/the-nnn-group/tutorials) were followed for US/WHAM, TAMD, TASS and bucket sampling.

## Acknowledgements

We thank Prof. Nisanth N. Nair for his guidance throughout this project.
