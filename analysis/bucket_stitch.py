#!/usr/bin/env python3
"""bucket_stitch.py - FES from 3 phi-restricted OPES buckets, reweighted and stitched.

Each bucket k: weights exp((opes.bias + wall)/kT) (OPES bias + the phi walls are the whole bias), histogram
on the common 72x72 grid -> F_k(phi,psi) up to a constant.  The buckets are windows on the phi circle that
overlap (b1/b2 on [-1.2,-0.8], b2/b3 on [0.8,1.2], b3/b1 on [2.94,3.34]); the 3 constants are fitted by least
squares on the overlap bins and the final F is the taper-weighted average of the shifted F_k.
t_eff = 3 x (time per bucket).   usage: bucket_stitch.py <dir_with_bucket_b1..b3> <prefix> <t_eff list ns>
"""
import sys, os
import numpy as np
import fesutils as fu

kT = fu.KB * 300.0
NB = fu.NBINS
_, C = fu.grid(NB)
# free arc of each bucket on the phi circle: centre, half-width (rad)
WIN = {1: (np.mean([2.94, -0.80 + 2*np.pi]), (-0.80 + 2*np.pi - 2.94) / 2),
       2: (0.0, 1.2),
       3: (np.mean([0.80, -2.94 + 2*np.pi]), (-2.94 + 2*np.pi - 0.80) / 2)}
TAPER = 0.4


def wmap(k):
    c, hw = WIN[k]
    d = np.abs(fu.wrap(C - c))
    return np.clip((hw - d) / TAPER, 0, 1), d < hw


def stitch(Fk):
    """Fk: {k: F (NBxNB, nan = unvisited)} -> stitched F (min 0)"""
    W = {k: wmap(k)[0][:, None] * np.ones(NB)[None, :] for k in Fk}
    ok = {k: np.isfinite(Fk[k]) & (Fk[k] < 40) for k in Fk}
    # least squares for offsets c_k (c_1 = 0): minimise sum w_i w_j (F_i + c_i - F_j - c_j)^2 over overlaps
    A, b = [], []
    for i in Fk:
        for j in Fk:
            if i >= j: continue
            m = ok[i] & ok[j] & (W[i] > 0) & (W[j] > 0)
            if not m.any(): continue
            ww = np.sqrt((W[i] * W[j])[m])
            d = (Fk[i] - Fk[j])[m]
            for wk, dk in zip(ww, d):
                row = np.zeros(3); row[i - 1] += wk; row[j - 1] -= wk
                A.append(row); b.append(-wk * dk)
    A = np.array(A); b = np.array(b)
    A, b = A[:, 1:], b                       # fix c_1 = 0
    c = np.r_[0.0, np.linalg.lstsq(A, b, rcond=None)[0]]
    num = np.zeros((NB, NB)); den = np.zeros((NB, NB))
    for k in Fk:
        v = np.isfinite(Fk[k]) & (W[k] > 0)
        num[v] += W[k][v] * (Fk[k][v] + c[k - 1]); den[v] += W[k][v]
    F = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
    return F - np.nanmin(F), c


def main():
    root, prefix, teffs = sys.argv[1], sys.argv[2], fu.parse_list(sys.argv[3])
    cv = {k: fu.load_colvar(os.path.join(root, f"bucket_b{k}", "COLVAR")) for k in (1, 2, 3)}
    for te in teffs:
        tb = te / 3.0 * 1000.0                         # ps per bucket
        Fk = {}
        for k in (1, 2, 3):
            s = fu.select_time(cv[k], 0.0, tb + 1e-6)
            w = fu.log_weights_to_weights((s["opes.bias"] + s["wall"]) / kT)
            H = fu.hist2d(s["phi"], s["psi"], NB, weights=w)
            # only keep bins inside the bucket window (outside it the wall pushes the walker away)
            _, inwin = wmap(k)
            H[~inwin, :] = 0
            F = fu.fes_from_prob(H / H.sum(), kT)
            Fk[k] = F
        F, c = stitch(Fk)
        out = f"{prefix}_{fu.fmt_ns(te)}ns.dat"
        fu.write_fes(out, F, header=f"stitched 3 buckets, {te/3:.2f} ns each, offsets {c.round(2)}")
        print(f"{out}: offsets(kJ/mol) {c.round(2)}  unvisited bins {int(np.isnan(F).sum())}")


if __name__ == "__main__":
    main()
