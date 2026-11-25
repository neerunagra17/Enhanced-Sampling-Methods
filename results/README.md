# Results (interim: four of five methods)

Free-energy surfaces and L² errors of umbrella sampling, TASS, WT-MetaD (two estimators) and
TAMD at the common t_eff checkpoints (10–500 ns). OPES results are pending.

Files
- `fes_<method>_<t>ns.dat` — FES on the common 72 × 72 grid, φ, ψ in rad, F in kJ/mol, min = 0,
  unvisited bins non-finite. `<method>`: `us`, `tass`, `metad` (bias-derived, sum of hills),
  `metadrw` (reweighted with `metad.rbias`), `tamd`.
- `l2_ref.csv`, `l2_external_table.csv` — L² against the GASS reference (mean alignment; the
  second file also gives minimum alignment).
- `l2_internal.csv` — L² against the same method's 500 ns surface.
- `l2_pairwise_200ns.csv`, `l2_pairwise_500ns.csv` — L² between methods; row = reference
  (defines the mask), column = estimate.
- `basin_features_500ns.csv` — basin free energies, barriers along φ and global minimum (kJ/mol).

Conventions (identical to the report): all surfaces shifted to min = 0; mask F_ref < 14.58
kcal/mol (61.0 kJ/mol); points the estimate has not visited are excluded and counted
(`missing`); F_sim mean-aligned to F_ref on the mask; L² in kcal/mol, surfaces in kJ/mol.
Convergence (`../analysis/convergence_table.py`): first checkpoint from which L² stays below
k_B T (0.596 kcal/mol) or 0.5 kcal/mol with at most 1 % of the mask unvisited at all later
checkpoints.

The bias-derived WT-MetaD surfaces are F = −Σ(stored hills): PLUMED writes well-tempered hill
heights already multiplied by γ/(γ−1), so no further factor is applied (see
`../analysis/sumhills_py.py`).

Regenerate:
    DATA=/path/to/production_output bash ../analysis/run_fes_checkpoints.sh
    PYTHONPATH=../analysis python3 ../analysis/plot_compare.py               # internal + pairwise
    PYTHONPATH=../analysis python3 ../analysis/plot_compare.py --ref fes_ref.dat
    PYTHONPATH=../analysis python3 ../analysis/plot_external.py               # needs fes_ref.dat
    python3 ../analysis/convergence_table.py l2_ref.csv
