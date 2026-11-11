#!/usr/bin/env python3
"""l2_error.py - L2 error of FES estimates with respect to a reference.

    L2 = sqrt( (1/N_g) * sum_i [F_ref(i) - F_sim(i)]^2 ),
    sum over grid points with F_ref(i) < F_cut   (N_g = number of such points)

Both surfaces are put on the common grid and shifted to min = 0 (or, with
--align mean, F_sim is shifted to have the same mean as F_ref over the
points used, which removes the arbitrary constant optimally).

Examples
  # error w.r.t. an independent reference
  python3 l2_error.py --ref ref/fes_ref.dat --sim fes_opes_*ns.dat \
          --fcut 40 --method OPES --out l2_opes_ref.csv
  # internal convergence: reference = the same method at its longest t_eff
  python3 l2_error.py --ref fes_opes_200ns.dat --sim fes_opes_*ns.dat \
          --fcut 40 --method OPES --out l2_opes_int.csv

t_eff is read from the file names (..._<t>ns.dat).
Use --sim-scale 4.184 if the simulated files are in kcal/mol.
"""
import argparse
import numpy as np
import fesutils as fu


def l2(Fref, Fsim, fcut, align="min", fill=None):
    mask = np.isfinite(Fref) & (Fref < fcut)
    s = Fsim.copy()
    n_missing = int(np.sum(mask & ~np.isfinite(s)))
    if n_missing:
        if fill is None:                      # ignore points the method never visited
            mask &= np.isfinite(s)
        else:                                 # or count them with a penalty value
            s[mask & ~np.isfinite(s)] = fill
    r, s = Fref[mask], s[mask]
    if align == "mean":
        s = s - s.mean() + r.mean()
    return np.sqrt(np.mean((r - s) ** 2)), int(mask.sum()), n_missing


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--sim", nargs="+", required=True)
    ap.add_argument("--fcut", type=float, required=True, help="kJ/mol")
    ap.add_argument("--align", choices=["min", "mean"], default="min")
    ap.add_argument("--fill", type=float, default=None,
                    help="value given to unvisited points inside the mask "
                         "(default: leave them out and report how many)")
    ap.add_argument("--ref-scale", type=float, default=1.0)
    ap.add_argument("--sim-scale", type=float, default=1.0)
    ap.add_argument("--nbins", type=int, default=fu.NBINS)
    ap.add_argument("--method", default="method")
    ap.add_argument("--out", required=True, help="CSV: method,t_eff_ns,L2,N_g,missing")
    a = ap.parse_args()

    Fref = fu.read_fes(a.ref, a.nbins, a.ref_scale)
    rows = []
    for path in a.sim:
        t = fu.teff_from_name(path)
        Fsim = fu.read_fes(path, a.nbins, a.sim_scale)
        e, ng, miss = l2(Fref, Fsim, a.fcut, a.align, a.fill)
        rows.append((t, e, ng, miss, path))
    rows.sort(key=lambda r: r[0])
    with open(a.out, "w") as fh:
        fh.write("method,t_eff_ns,L2_kJmol,N_g,missing,file\n")
        for t, e, ng, miss, path in rows:
            fh.write(f"{a.method},{t:g},{e:.4f},{ng},{miss},{path}\n")
            print(f"{a.method:10s} t_eff={t:8g} ns  L2={e:7.3f} kJ/mol  "
                  f"N_g={ng}  unvisited={miss}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
