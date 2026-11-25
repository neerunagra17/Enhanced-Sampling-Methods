#!/usr/bin/env python3
"""plot_traces.py - CV time series; used for two checks.

1. Did the run cross the phi barrier? (MetaD, OPES, OPES+bucket)
     python3 plot_traces.py --colvar COLVAR --cols phi,psi --out traces_opes.png
2. TAMD/TASS adiabatic check: does the physical CV q follow the auxiliary s?
     python3 plot_traces.py --colvar COLVAR_cc --cols phi \
             --colvar2 COLVAR --cols2 ext.phi_fict --tmax 200 --out q_vs_s.png
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fesutils as fu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--colvar", required=True)
    ap.add_argument("--cols", required=True, help="comma list of columns")
    ap.add_argument("--colvar2", default=None, help="second file (overlay)")
    ap.add_argument("--cols2", default=None)
    ap.add_argument("--tmax", type=float, default=None, help="ps (zoom)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cols = a.cols.split(",")
    fig, axes = plt.subplots(len(cols), 1, figsize=(8, 2.4 * len(cols)), squeeze=False)
    cv = fu.select_time(fu.load_colvar(a.colvar), None, a.tmax)
    cv2 = None
    if a.colvar2:
        cv2 = fu.select_time(fu.load_colvar(a.colvar2), None, a.tmax)
        cols2 = a.cols2.split(",")
    for k, c in enumerate(cols):
        ax = axes[k, 0]
        ax.plot(cv["time"] / 1000, fu.wrap(cv[c]), lw=0.5, label=c)
        if cv2 is not None:
            ax.plot(cv2["time"] / 1000, fu.wrap(cv2[cols2[k]]), lw=0.5,
                    alpha=0.8, label=cols2[k])
            ax.legend(frameon=False, loc="upper right", fontsize=8)
        ax.set_ylabel(f"{c} (rad)")
    axes[-1, 0].set_xlabel("time (ns)")
    fig.tight_layout()
    fig.savefig(a.out, dpi=200)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
