# Reference surface

The external reference used in the report is a free-energy surface F(φ, ψ) of alanine
dipeptide in water computed with globally accelerated sliced sampling (GASS). It was
provided for this project and is **not redistributed** here.

To reproduce the external comparison, place the file in `results/fes_ref.dat`
(three columns `phi psi F`, angles in rad, **F in kcal/mol**; non-finite entries such as
`Infinity` are treated as undefined). `plot_compare.py --ref` and `plot_external.py` read it
with `fesutils.read_fes(..., energy_scale=4.184)` and interpolate it periodically
(bilinear) onto the common 72 × 72 grid.

The simulation conditions of the GASS calculation are not recorded with the file, so the
surface is used as an external benchmark, not as an exact reference (report, Section 4.4).
