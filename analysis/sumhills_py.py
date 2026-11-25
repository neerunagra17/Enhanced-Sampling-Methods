#!/usr/bin/env python3
"""sumhills_py.py - pure-python replacement for 'plumed sum_hills' (2D WT-MetaD).

Sums the Gaussians in HILLS deposited before each t_eff and writes the free
energy on the common 72x72 grid (periodic, min shifted to 0), in the common FES
format. Used because PLUMED is not installed on the analysis machine.

For well-tempered MetaD, PLUMED (MetaD.cpp, 2.9.2) writes each hill height in
HILLS already multiplied by gamma/(gamma-1). The sum of the stored Gaussians is
therefore -F directly, exactly as 'plumed sum_hills' returns it:
    F = -sum(stored Gaussians)
No further gamma/(gamma-1) factor may be applied.
    python3 sumhills_py.py --hills HILLS --teff 10,20,50 --prefix fes_metad
"""
import argparse
import numpy as np
import fesutils as fu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hills", required=True)
    ap.add_argument("--teff", required=True)
    ap.add_argument("--prefix", default="fes_metad")
    ap.add_argument("--nbins", type=int, default=fu.NBINS)
    a = ap.parse_args()
    d = np.loadtxt(a.hills, comments="#")
    t = fu.continuous_time(d[:, 0])
    phi, psi, sp, ss, h, gam = d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5], d[:, 6]
    print(f"{len(t)} hills, last at {t[-1]/1000:.2f} ns, biasfactor {gam[0]:g}")
    _, c = fu.grid(a.nbins)

    def g(x0, s):                              # (N, nbins) periodic gaussian
        dx = fu.wrap(c[None, :] - x0[:, None])
        return np.exp(-0.5 * (dx / s[:, None]) ** 2)

    for t_ns in fu.parse_list(a.teff):
        n = int(np.searchsorted(t, t_ns * 1000.0, side="right"))
        V = np.zeros((a.nbins, a.nbins))
        for i in range(0, n, 20000):
            s = slice(i, min(i + 20000, n))
            V += (g(phi[s], sp[s]) * h[s, None]).T @ g(psi[s], ss[s])
        F = -V                                 # heights already include gamma/(gamma-1)
        F -= F.min()
        out = f"{a.prefix}_{fu.fmt_ns(t_ns)}ns.dat"
        fu.write_fes(out, F, header=f"sum_hills(py) {n} hills")
        print(f"  wrote {out} ({n} hills)")


if __name__ == "__main__":
    main()
