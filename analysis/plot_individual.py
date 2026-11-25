#!/usr/bin/env python3
"""plot_individual.py - one set of files per method: 2D FES, 1D profiles, L2 convergence.
Run inside results/ (PYTHONPATH=../analysis). Output: individual/<method>/...
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plot_compare as pc
from plot_compare import files, load, proj, C, CMAP, COL, METHODS, FCUT, l2k, KCAL
import fesutils as fu

for m, name in METHODS:
    fs = files(m)
    d = f"individual/{m}"; os.makedirs(d, exist_ok=True)
    tref = max(fs)
    for t in (200.0, 500.0):
        if t not in fs: continue
        F = load(fs[t])
        # 2D
        fig, ax = plt.subplots(figsize=(5.4, 4.5))
        lev = np.arange(0, 65, 5)
        cf = ax.contourf(C, C, np.clip(F, 0, 60).T, levels=lev, cmap=CMAP, extend="max")
        ax.contour(C, C, np.clip(F, 0, 60).T, levels=lev, colors="k", linewidths=0.25)
        fig.colorbar(cf, label=r"$F$ (kJ/mol)")
        ax.set_xlabel(r"$\phi$ (rad)"); ax.set_ylabel(r"$\psi$ (rad)"); ax.set_aspect("equal")
        ax.set_title(f"{name}, {t:g} ns"); fig.tight_layout()
        fig.savefig(f"{d}/fes2d_{t:g}ns.png", dpi=170); fig.savefig(f"{d}/fes2d_{t:g}ns.pdf"); plt.close(fig)
        # 1D
        fig, axs = plt.subplots(1, 2, figsize=(9, 3.5))
        for ax, ai, lab in zip(axs, (1, 0), (r"$\phi$", r"$\psi$")):
            ax.plot(C, proj(F, ai), color=CLAY_LINE if (CLAY_LINE := "#5d7fb0") else None, lw=2)
            ax.set_xlabel(lab + " (rad)"); ax.set_ylabel(r"$F$ (kJ/mol)"); ax.grid(alpha=.3)
            ax.set_xlim(-np.pi, np.pi); ax.set_ylim(0, 65)
        fig.suptitle(f"{name}: 1D profiles, {t:g} ns"); fig.tight_layout()
        fig.savefig(f"{d}/fes1d_{t:g}ns.png", dpi=170); fig.savefig(f"{d}/fes1d_{t:g}ns.pdf"); plt.close(fig)
    # L2 convergence vs own 500 ns
    Fr = load(fs[tref]); ts, es = [], []
    for t, p in fs.items():
        if t >= tref: continue
        ts.append(t); es.append(l2k(Fr, load(p))[0])
    fig, ax = plt.subplots(figsize=(5.4, 4))
    ax.plot(ts, es, "o-", color="#c5559b", lw=2, ms=5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(alpha=.3, which="both")
    ax.set_xlabel(r"$t_{\rm eff}$ (ns)"); ax.set_ylabel(r"$L^2$ error (kcal/mol)")
    ax.set_title(f"{name}: convergence vs own {tref:g} ns"); fig.tight_layout()
    fig.savefig(f"{d}/l2_internal.png", dpi=170); fig.savefig(f"{d}/l2_internal.pdf"); plt.close(fig)
    np.savetxt(f"{d}/l2_internal.csv", np.c_[ts, es], delimiter=",", header="t_eff_ns,L2_kcalmol", fmt="%.4f")
    if os.path.exists("fes_ref.dat"):
        Fx = load("fes_ref.dat", KCAL); tx, ex = [], []
        for t, p in fs.items():
            tx.append(t); ex.append(l2k(Fx, load(p))[0])
        fig, ax = plt.subplots(figsize=(5.4, 4))
        ax.plot(tx, ex, "o-", color="#5d7fb0", lw=2, ms=5)
        ax.set_xscale("log"); ax.grid(alpha=.3, which="both")
        ax.set_xlabel(r"$t_{\rm eff}$ (ns)"); ax.set_ylabel(r"$L^2$ error (kcal/mol)")
        ax.set_title(f"{name}: error vs GASS reference"); fig.tight_layout()
        fig.savefig(f"{d}/l2_external.png", dpi=170); fig.savefig(f"{d}/l2_external.pdf"); plt.close(fig)
        np.savetxt(f"{d}/l2_external.csv", np.c_[tx, ex], delimiter=",", header="t_eff_ns,L2_kcalmol", fmt="%.4f")
    print(m, "done")
