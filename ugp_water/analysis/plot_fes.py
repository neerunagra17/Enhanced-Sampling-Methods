#!/usr/bin/env python3
"""plot_fes.py - 2D FES figure (common format) for the report.

    python3 plot_fes.py fes_opes_200ns.dat --title "OPES, 200 ns" --out fes_opes.pdf
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fesutils as fu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fes")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--fmax", type=float, default=60.0, help="colour cap, kJ/mol")
    ap.add_argument("--step", type=float, default=5.0, help="contour spacing, kJ/mol")
    ap.add_argument("--scale", type=float, default=1.0, help="4.184 for kcal/mol files")
    ap.add_argument("--nbins", type=int, default=fu.NBINS)
    a = ap.parse_args()

    F = fu.read_fes(a.fes, a.nbins, a.scale)
    _, c = fu.grid(a.nbins)
    Fp = np.clip(F, 0, a.fmax)                       # nan (unvisited) stays white
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    levels = np.arange(0, a.fmax + a.step, a.step)
    cf = ax.contourf(c, c, Fp.T, levels=levels, cmap="viridis")
    ax.contour(c, c, Fp.T, levels=levels, colors="k", linewidths=0.3)
    fig.colorbar(cf, ax=ax, label=r"$F$ (kJ mol$^{-1}$)")
    ax.set_xlabel(r"$\phi$ (rad)")
    ax.set_ylabel(r"$\psi$ (rad)")
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-np.pi, np.pi)
    ax.set_aspect("equal")
    if a.title:
        ax.set_title(a.title)
    fig.tight_layout()
    fig.savefig(a.out)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
