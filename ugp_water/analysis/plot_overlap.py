#!/usr/bin/env python3
"""plot_overlap.py - umbrella / TASS window diagnostics.

Top: histogram of the restrained CV for every window (overlap check).
Bottom: the restrained CV vs time for every window (equilibration check).

    python3 plot_overlap.py --meta meta.dat --tmin 500 --out overlap.pdf
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fesutils as fu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--phi", default="phi")
    ap.add_argument("--tmin", type=float, default=0.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 6.5))
    edges = np.linspace(-np.pi, np.pi, 361)
    centres = 0.5 * (edges[1:] + edges[:-1])
    paths = [l.split()[0] for l in open(a.meta) if l.strip() and not l.startswith("#")]
    cmap = plt.cm.viridis(np.linspace(0, 1, len(paths)))
    for col, path in zip(cmap, paths):
        cv = fu.select_time(fu.load_colvar(path), a.tmin, None)
        h, _ = np.histogram(fu.wrap(cv[a.phi]), bins=edges, density=True)
        ax1.plot(centres, h, color=col, lw=1)
        ax2.plot(cv["time"] / 1000, fu.wrap(cv[a.phi]), ",", color=col)
    ax1.set_xlabel(r"$\phi$ (rad)")
    ax1.set_ylabel("probability density")
    ax1.set_xlim(-np.pi, np.pi)
    ax2.set_xlabel("time (ns)")
    ax2.set_ylabel(r"$\phi$ (rad)")
    fig.tight_layout()
    fig.savefig(a.out)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
