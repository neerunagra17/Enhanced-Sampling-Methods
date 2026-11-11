#!/usr/bin/env python3
"""plot_l2.py - the main result: L2 error vs t_eff, one curve per method.

    python3 plot_l2.py --csv l2_us_ref.csv l2_metad_ref.csv l2_opes_ref.csv \
            --title "Error w.r.t. reference" --out l2_vs_teff_ref.pdf
"""
import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", nargs="+", required=True)
    ap.add_argument("--out", required=True, help=".pdf (vector, for the report) or .png")
    ap.add_argument("--title", default="")
    ap.add_argument("--logy", action="store_true")
    a = ap.parse_args()

    markers = "osD^vP*X"
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for k, path in enumerate(a.csv):
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        t = [float(r["t_eff_ns"]) for r in rows]
        e = [float(r["L2_kJmol"]) for r in rows]
        ax.plot(t, e, marker=markers[k % len(markers)], lw=1.6, ms=5,
                label=rows[0]["method"])
    ax.set_xlabel(r"$t_{\mathrm{eff}}$ (ns)")
    ax.set_ylabel(r"$L^2$ error (kJ mol$^{-1}$)")
    if a.logy:
        ax.set_yscale("log")
    if a.title:
        ax.set_title(a.title)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(a.out)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
