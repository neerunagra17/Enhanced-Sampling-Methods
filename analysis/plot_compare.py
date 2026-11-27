#!/usr/bin/env python3
"""plot_compare.py - 2D FES, 1D projections and L2 comparison for all methods.

Run inside results/:   python3 ../analysis/plot_compare.py [--ref fes_ref.dat]
Without --ref: internal convergence (each method vs its own longest run) and the
pairwise L2 matrix between methods. With --ref: L2 vs t_eff against the reference.
"""
import argparse, glob, re, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fesutils as fu
from l2_error import l2

KT = fu.KB * 300.0
from matplotlib.colors import LinearSegmentedColormap
RAINBOW = [(0.00, "#ff0000"), (0.10, "#ff9a00"), (0.20, "#ffee00"), (0.30, "#b6f000"), (0.40, "#10e010"),
           (0.52, "#00a8a0"), (0.60, "#0a28ff"), (0.72, "#a020f0"), (0.82, "#ff10e0"), (1.00, "#ff1070")]  # red -> pink
CLAY = [c for _, c in RAINBOW]
CMAP = LinearSegmentedColormap.from_list("rainbow_red_pink", RAINBOW)
METHODS = [("us", "Umbrella sampling"), ("tass", "TASS"), ("metad", "WT-MetaD (sum hills)"),
           ("metadrw", "WT-MetaD (reweighted)"), ("tamd", "TAMD"),
           ("opes", "OPES"), ("opesb", "OPES + bucket")]
COL = {"us": "C0", "tass": "C1", "metad": "C2", "metadrw": "C3", "tamd": "C4", "opes": "C5", "opesb": "C6"}
KCAL = 4.184
FCUT, NB = 14.58 * KCAL, fu.NBINS     # 14.58 kcal/mol cutoff (= 61.0 kJ/mol); L2 reported in kcal/mol

_, C = fu.grid(NB)


def files(m):
    out = {}
    for p in glob.glob(f"fes_{m}_*ns.dat"):
        mm = re.fullmatch(rf"fes_{m}_([0-9.]+)ns\.dat", p)
        if mm:
            out[float(mm.group(1))] = p
    return dict(sorted(out.items()))


def load(p, scale=1.0):
    return fu.read_fes(p, NB, scale)


def l2k(Fr, Fs):
    """L2 in kcal/mol: mask = reference < 14.58 kcal/mol and simulation visited; mean-aligned."""
    e, ng, miss = l2(Fr, Fs, FCUT, "mean")
    return e / KCAL, ng, miss


def proj(F, axis):
    """1D free energy: integrate Boltzmann factor over the other CV (nan = never visited)."""
    P = np.where(np.isfinite(F), np.exp(-F / KT), 0.0).sum(axis=axis)
    with np.errstate(divide="ignore"):
        f = -KT * np.log(P)
    return f - np.nanmin(f)


def fes2d(tag, t, ref_label=None):
    fig, axs = plt.subplots(1, len(METHODS), figsize=(3.6 * len(METHODS), 3.5), sharey=True)
    lev = np.arange(0, 65, 5)
    for ax, (m, name) in zip(axs, METHODS):
        fs = files(m)
        if t not in fs:
            ax.set_title(f"{name}\n(no {t:g} ns)"); continue
        F = np.clip(load(fs[t]), 0, 60)
        cf = ax.contourf(C, C, F.T, levels=lev, cmap=CMAP, extend="max")
        ax.contour(C, C, F.T, levels=lev, colors="k", linewidths=0.25)
        ax.set_title(f"{name}, {t:g} ns", fontsize=9); ax.set_xlabel(r"$\phi$ (rad)")
        ax.set_aspect("equal")
    axs[0].set_ylabel(r"$\psi$ (rad)")
    fig.colorbar(cf, ax=axs, label=r"$F$ (kJ/mol)", shrink=0.85, pad=0.01)
    fig.savefig(f"fes2d_{tag}.pdf", bbox_inches="tight"); fig.savefig(f"fes2d_{tag}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fes1d(tag, t):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    for m, name in METHODS:
        fs = files(m)
        if t not in fs: continue
        F = load(fs[t])
        for ax, ax_i in zip(axs, (1, 0)):
            ax.plot(C, proj(F, ax_i), color=COL[m], lw=1.6, label=name)
    for ax, lab in zip(axs, (r"$\phi$", r"$\psi$")):
        ax.set_xlabel(lab + " (rad)"); ax.set_ylabel(r"$F$ (kJ/mol)"); ax.grid(alpha=.3)
        ax.set_xlim(-np.pi, np.pi); ax.set_ylim(0, 65)
    axs[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"1D free energy profiles (projected), {t:g} ns")
    fig.tight_layout(); fig.savefig(f"fes1d_{tag}.pdf"); fig.savefig(f"fes1d_{tag}.png", dpi=150)
    plt.close(fig)


def convergence(tref, ref_path=None, tag="internal", title=""):
    fig, ax = plt.subplots(figsize=(6.3, 4.4))
    rows = []
    for m, name in METHODS:
        fs = files(m)
        if not fs: continue
        Fr = load(ref_path, KCAL) if ref_path else (load(fs[tref]) if tref in fs else None)
        if Fr is None: continue
        ts, es = [], []
        for t, p in fs.items():
            if not ref_path and t >= tref: continue       # skip the trivial 0 at the reference itself
            e, ng, miss = l2k(Fr, load(p))
            ts.append(t); es.append(e); rows.append((name, t, e, ng, miss))
        ax.plot(ts, es, "o-", color=COL[m], lw=1.6, ms=4, label=name)
    ax.set_xlabel(r"$t_{\rm eff}$ (ns)"); ax.set_ylabel(r"$L^2$ error (kcal/mol)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(alpha=.3, which="both"); ax.legend(frameon=False, fontsize=8)
    ax.set_title(title, fontsize=10); fig.tight_layout()
    fig.savefig(f"l2_{tag}.pdf"); fig.savefig(f"l2_{tag}.png", dpi=150); plt.close(fig)
    with open(f"l2_{tag}.csv", "w") as fh:
        fh.write("method,t_eff_ns,L2_kcalmol,N_g,missing\n")
        for r in rows: fh.write(f"{r[0]},{r[1]:g},{r[2]:.4f},{r[3]},{r[4]}\n")


def pairwise(t):
    ms = [(m, n) for m, n in METHODS if t in files(m)]
    Fs = [load(files(m)[t]) for m, _ in ms]
    M = np.array([[l2k(a, b)[0] for b in Fs] for a in Fs])   # row = reference
    fig, ax = plt.subplots(figsize=(5.8, 4.8))
    im = ax.imshow(M, cmap="magma_r")
    ax.set_xticks(range(len(ms))); ax.set_yticks(range(len(ms)))
    ax.set_xticklabels([n for _, n in ms], rotation=35, ha="right", fontsize=8)
    ax.set_yticklabels([n for _, n in ms], fontsize=8)
    for i in range(len(ms)):
        for j in range(len(ms)):
            ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=8,
                    color="w" if M[i, j] > M.max() / 2 else "k")
    fig.colorbar(im, label=r"$L^2$ (kcal/mol)"); ax.set_title(f"Pairwise $L^2$ (row = reference), {t:g} ns")
    fig.tight_layout(); fig.savefig(f"l2_pairwise_{t:g}ns.pdf"); fig.savefig(f"l2_pairwise_{t:g}ns.png", dpi=150)
    plt.close(fig)
    np.savetxt(f"l2_pairwise_{t:g}ns.csv", M, delimiter=",", header=",".join(n for _, n in ms), fmt="%.4f")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--ref", default=None)
    a = ap.parse_args()
    for t in (200.0, 500.0):
        fes2d(f"{t:g}ns", t); fes1d(f"{t:g}ns", t); pairwise(t)
    convergence(500.0, None, "internal", "Internal convergence (vs each method's 500 ns)")
    if a.ref:
        convergence(None, a.ref, "ref", "Error w.r.t. reference")
