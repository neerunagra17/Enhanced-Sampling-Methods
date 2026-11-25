#!/usr/bin/env python3
"""plot_external.py - comparison with the external (GASS) reference FES, in kcal/mol.
Run inside results/ (PYTHONPATH=../analysis); needs fes_ref.dat (kcal/mol, any grid).
Reference is resampled (periodic bilinear) onto our 72x72 grid, both surfaces shifted to
min=0; L2 over points with F_ref < 14.58 kcal/mol that the simulation visited, mean-aligned.
Outputs: ext_fes2d_500ns, ext_diff_500ns, ext_fes1d_500ns, l2_external_table.csv
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_compare import files, load, proj, C, CMAP, COL, METHODS, FCUT, KCAL, l2k, KT
from l2_error import l2

R = load("fes_ref.dat", KCAL)           # kJ/mol on our grid
T = 500.0
ms = [(m, n) for m, n in METHODS if T in files(m)]
Fs = {m: load(files(m)[T]) for m, _ in ms}
lev = np.arange(0, 65, 5)

fig, axs = plt.subplots(1, len(ms) + 1, figsize=(3.4 * (len(ms) + 1), 3.5), sharey=True)
for ax, (nm, F) in zip(axs, [("GASS reference", R)] + [(n, Fs[m]) for m, n in ms]):
    cf = ax.contourf(C, C, np.clip(F, 0, 60).T, levels=lev, cmap=CMAP, extend="max")
    ax.contour(C, C, np.clip(F, 0, 60).T, levels=lev, colors="k", linewidths=0.25)
    ax.set_title(nm if nm.startswith("GASS") else f"{nm}, {T:g} ns", fontsize=9)
    ax.set_xlabel(r"$\phi$ (rad)"); ax.set_aspect("equal")
axs[0].set_ylabel(r"$\psi$ (rad)")
fig.colorbar(cf, ax=axs, label=r"$F$ (kJ/mol)", shrink=0.85, pad=0.01)
fig.savefig("ext_fes2d_500ns.png", dpi=150, bbox_inches="tight"); fig.savefig("ext_fes2d_500ns.pdf", bbox_inches="tight"); plt.close(fig)

# difference maps (mean-aligned over the L2 mask), kcal/mol
fig, axs = plt.subplots(1, len(ms), figsize=(3.6 * len(ms), 3.5), sharey=True)
for ax, (m, n) in zip(np.atleast_1d(axs), ms):
    F = Fs[m]; mask = np.isfinite(R) & (R < FCUT) & np.isfinite(F)
    D = np.full(R.shape, np.nan); off = (F[mask] - R[mask]).mean()
    D[mask] = (F[mask] - off - R[mask]) / KCAL
    im = ax.pcolormesh(C, C, D.T, cmap="RdBu_r", vmin=-3, vmax=3, shading="nearest")
    ax.set_title(f"{n}\nL2={l2k(R, F)[0]:.2f} kcal/mol", fontsize=9); ax.set_xlabel(r"$\phi$ (rad)"); ax.set_aspect("equal")
np.atleast_1d(axs)[0].set_ylabel(r"$\psi$ (rad)")
fig.colorbar(im, ax=axs, label="F_sim - F_ref (kcal/mol)", shrink=0.85, pad=0.01)
fig.savefig("ext_diff_500ns.png", dpi=150, bbox_inches="tight"); fig.savefig("ext_diff_500ns.pdf", bbox_inches="tight"); plt.close(fig)

# 1D profiles vs reference
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
for ax, ai, lab in zip(axs, (1, 0), (r"$\phi$", r"$\psi$")):
    ax.plot(C, proj(R, ai), "k--", lw=2, label="GASS reference")
    for m, n in ms: ax.plot(C, proj(Fs[m], ai), color=COL[m], lw=1.4, label=n)
    ax.set_xlabel(lab + " (rad)"); ax.set_ylabel(r"$F$ (kJ/mol)"); ax.grid(alpha=.3)
    ax.set_xlim(-np.pi, np.pi); ax.set_ylim(0, 65)
axs[0].legend(frameon=False, fontsize=8)
fig.suptitle(f"1D profiles vs reference, {T:g} ns"); fig.tight_layout()
fig.savefig("ext_fes1d_500ns.png", dpi=150); fig.savefig("ext_fes1d_500ns.pdf"); plt.close(fig)

# table: both alignments, several t_eff
with open("l2_external_table.csv", "w") as fh:
    fh.write("method,t_eff_ns,L2_mean_aligned_kcal,L2_min_aligned_kcal,N_g,missing\n")
    for m, n in METHODS:
        for t, p in files(m).items():
            F = load(p); e1, ng, miss = l2k(R, F)
            e2 = l2(R, F, FCUT, "min")[0] / KCAL
            fh.write(f"{n},{t:g},{e1:.4f},{e2:.4f},{ng},{miss}\n")
print(open("l2_external_table.csv").read())
