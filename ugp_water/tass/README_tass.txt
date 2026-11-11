TASS = umbrella windows (like Section 2) + extended variables at T_s + WT-MetaD on psi.

1. Preparation: identical to umbrella sampling. Reuse Section 2 scripts:
      (windows.sh, us_prep.in and plumed_us.template are already in tass/)
      bash ../us/run_us_prep.sh          (physical phi restrained, 100 ps per window)
2. Production (window = folder win_<centre>):
      TEMPLATE=plumed_tass.template MDIN=tass_prod.in bash ../us/run_us_prod.sh
3. Reconstruction:
   a) official: the TASS reweighting / WHAM (or mean-force, TASS-MF) code from
      Prof. Nair's GitHub - follow its README.
   b) cross-check with the supplied wham2d.py (weighted frames, WHAM at T_s):
      python3 ../analysis/wham2d.py --meta meta.dat --phi ext.phi_fict \
          --psi ext.psi_fict --logw metad.rbias --T 900 --tmin 500 --out fes_tass.dat
t_eff = (number of windows) x (time per window), as for umbrella sampling.
