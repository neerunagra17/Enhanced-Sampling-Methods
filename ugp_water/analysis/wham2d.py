#!/usr/bin/env python3
"""wham2d.py - 2D WHAM: windows restrained along phi only, FES in (phi, psi).

Cross-check for (and fallback to) the WHAM code supplied by the group.
Also used for TASS (weighted frames, auxiliary temperature); see Section 6.

Metafile: one line per window
    path/to/COLVAR   centre(rad)   kappa(kJ/mol/rad^2)
Lines starting with '#' are ignored.

Examples
  # all data after 500 ps of equilibration in every window
  python3 wham2d.py --meta meta.dat --tmin 500 --out fes_us.dat
  # only the first 2 ns of production per window (for t_eff checkpoints)
  python3 wham2d.py --meta meta.dat --tmin 500 --tw 2000 --out fes_us_64ns.dat

Theory: Eqs. (WHAM) and (f_i) of the report. The bias depends on phi only,
so the WHAM equations are solved on a fine phi grid (--nphi-fine) to avoid
discretisation error, and P(phi, psi) is then summed onto the common grid.
"""
import argparse
import numpy as np
from scipy.special import logsumexp
import fesutils as fu


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta", required=True, help="metafile (see top of file)")
    ap.add_argument("--out", required=True, help="output FES (common format)")
    ap.add_argument("--T", type=float, default=300.0,
                    help="temperature (K) used in WHAM and F=-kT lnP "
                         "(TASS: the auxiliary temperature)")
    ap.add_argument("--phi", default="phi", help="COLVAR column of the restrained CV")
    ap.add_argument("--psi", default="psi", help="COLVAR column of the free CV")
    ap.add_argument("--logw", default=None,
                    help="optional column with a per-frame reweighting bias "
                         "(kJ/mol), frame weight = exp(+bias/kT) (TASS: metad.rbias)")
    ap.add_argument("--tmin", type=float, default=0.0, help="discard t < tmin (ps)")
    ap.add_argument("--tw", type=float, default=None,
                    help="use only tmin <= t < tmin+tw (ps) of every window")
    ap.add_argument("--nbins", type=int, default=fu.NBINS, help="output grid")
    ap.add_argument("--nphi-fine", type=int, default=720,
                    help="phi bins used internally (multiple of --nbins)")
    ap.add_argument("--tol", type=float, default=1e-6, help="convergence (kJ/mol)")
    ap.add_argument("--maxiter", type=int, default=100000)
    a = ap.parse_args()

    if a.nphi_fine % a.nbins:
        raise SystemExit("--nphi-fine must be a multiple of --nbins")
    kT = fu.KB * a.T
    beta = 1.0 / kT
    fine_edges = np.linspace(-np.pi, np.pi, a.nphi_fine + 1)
    fine_c = 0.5 * (fine_edges[1:] + fine_edges[:-1])
    psi_edges, _ = fu.grid(a.nbins)

    # ---------------- read all windows ----------------
    H, logN, W = [], [], []
    for line in open(a.meta):
        s = line.split()
        if not s or s[0].startswith("#"):
            continue
        path, centre, kappa = s[0], float(s[1]), float(s[2])
        cv = fu.load_colvar(path)
        tmax = None if a.tw is None else a.tmin + a.tw
        if tmax is not None and cv["time"][-1] < tmax - 1.0:
            print(f"  WARNING {path}: only {cv['time'][-1]:.0f} ps simulated, "
                  f"{tmax:.0f} ps requested - t_eff of this checkpoint is overstated")
        cv = fu.select_time(cv, a.tmin, tmax)
        n = len(cv[a.phi])
        if n == 0:
            raise SystemExit(f"{path}: no frames left after time selection")
        w = None
        if a.logw:
            w = fu.log_weights_to_weights(beta * cv[a.logw])
            w *= n / w.sum()                     # weighted counts sum to n
        h, _, _ = np.histogram2d(fu.wrap(cv[a.phi]), fu.wrap(cv[a.psi]),
                                 bins=[fine_edges, psi_edges], weights=w)
        d = fu.wrap(fine_c - centre)             # minimum-image distance
        H.append(h)
        logN.append(np.log(n))
        W.append(0.5 * kappa * d ** 2)
        print(f"  {path:40s} centre {centre:+7.3f}  frames {n}")
    H, logN, W = np.array(H), np.array(logN), np.array(W)   # W: (nwin, nphi)
    nwin = len(H)
    Htot = H.sum(axis=0)                                    # (nphi, npsi)
    print(f"{nwin} windows, {int(Htot.sum())} frames in total")

    # ---------------- self-consistent WHAM iterations ----------------
    f = np.zeros(nwin)
    with np.errstate(divide="ignore"):
        logHtot = np.log(Htot)
    for it in range(1, a.maxiter + 1):
        # log of sum_i N_i exp(beta (f_i - W_i(phi)))   -> shape (nphi,)
        logD = logsumexp(logN[:, None] + beta * (f[:, None] - W), axis=0)
        logP = logHtot - logD[:, None]                       # (nphi, npsi)
        # f_i = -kT ln sum_{phi,psi} P exp(-beta W_i)
        logPphi = logsumexp(logP, axis=1)                    # sum over psi
        fnew = -kT * logsumexp(logPphi[None, :] - beta * W, axis=1)
        fnew -= fnew[0]
        delta = np.max(np.abs(fnew - f))
        f = fnew
        if delta < a.tol:
            break
    print(f"WHAM converged in {it} iterations (max change {delta:.2e} kJ/mol)")

    # ---------------- sum fine phi bins onto the output grid ----------------
    P = np.exp(logP - np.nanmax(logP[np.isfinite(logP)]))
    P[~np.isfinite(P)] = 0.0
    P = P.reshape(a.nbins, a.nphi_fine // a.nbins, a.nbins).sum(axis=1)
    F = fu.fes_from_prob(P / P.sum(), kT)
    fu.write_fes(a.out, F, header=f"WHAM {nwin} windows T={a.T}K")
    np.savetxt(a.out + ".fwin", f, header="window free-energy offsets f_i (kJ/mol)")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
