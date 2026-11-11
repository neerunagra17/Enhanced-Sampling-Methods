#!/usr/bin/env python3
"""reweight_fes.py - FES from a single biased trajectory, at several t_eff.

One script for every single-trajectory method:

  OPES        weights exp(+V/kT), T = 300 K
      python3 reweight_fes.py --colvar COLVAR --bias opes.bias --T 300 \
              --teff 10,20,50,100 --prefix fes_opes
  MetaD (alternative to sum_hills; needs CALC_RCT in plumed.dat)
      python3 reweight_fes.py --colvar COLVAR --bias metad.rbias --T 300 \
              --teff 10,20,50,100 --prefix fes_metadrw --tmin 1000
  TAMD        no bias, histogram of the auxiliary variables at T_s
      python3 reweight_fes.py --colvar COLVAR --phi ext.phi_fict \
              --psi ext.psi_fict --T 1500 --teff 10,20,50 --prefix fes_tamd

For each value t in --teff (ns) it uses frames with tmin <= time < t
and writes  <prefix>_<t>ns.dat  in the common FES format (kJ/mol).
"""
import argparse
import numpy as np
import fesutils as fu


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--colvar", required=True)
    ap.add_argument("--phi", default="phi", help="COLVAR column for phi")
    ap.add_argument("--psi", default="psi", help="COLVAR column for psi")
    ap.add_argument("--bias", default=None,
                    help="COLVAR column with the bias (kJ/mol); omit for TAMD")
    ap.add_argument("--T", type=float, default=300.0,
                    help="temperature in kT (TAMD: auxiliary temperature T_s)")
    ap.add_argument("--teff", required=True, help="comma list of times in ns")
    ap.add_argument("--tmin", type=float, default=0.0,
                    help="discard frames before tmin (ps)")
    ap.add_argument("--nbins", type=int, default=fu.NBINS)
    ap.add_argument("--prefix", required=True)
    a = ap.parse_args()

    kT = fu.KB * a.T
    cv = fu.load_colvar(a.colvar)
    tlast = cv["time"][-1]
    print(f"{a.colvar}: {len(cv['time'])} frames, last time {tlast/1000:.2f} ns")

    for t_ns in fu.parse_list(a.teff):
        dt = cv["time"][1] - cv["time"][0] if len(cv["time"]) > 1 else 0.0
        if t_ns * 1000.0 > tlast + dt + 1e-6:
            print(f"  skip {t_ns:g} ns (trajectory too short)")
            continue
        sub = fu.select_time(cv, a.tmin, t_ns * 1000.0 + 1e-6)
        w = None
        if a.bias:
            w = fu.log_weights_to_weights(sub[a.bias] / kT)
        H = fu.hist2d(sub[a.phi], sub[a.psi], a.nbins, weights=w)
        F = fu.fes_from_prob(H / H.sum(), kT)
        out = f"{a.prefix}_{fu.fmt_ns(t_ns)}ns.dat"
        fu.write_fes(out, F, header=f"frames={len(sub['time'])} T={a.T}K")
        neff = "" if w is None else f"  effective samples {w.sum()**2/(w**2).sum():.0f}"
        print(f"  wrote {out}{neff}")


if __name__ == "__main__":
    main()
