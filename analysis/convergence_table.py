#!/usr/bin/env python3
"""convergence_table.py - first t_eff at which an L2 curve has converged.

Criterion (Tables 6.3 and 6.4 of the report): the first checkpoint from which L2
stays below the threshold AND at most 1 % of the mask points are unvisited at that
checkpoint and at every later one. Thresholds: k_B T at 300 K (0.596 kcal/mol)
and 0.5 kcal/mol. "--" means not reached within the run.

Input: a CSV written by plot_compare.py (l2_ref.csv or l2_internal.csv), columns
method,t_eff_ns,L2_kcalmol,N_g,missing. The mask size is N_g + missing.

    python3 convergence_table.py l2_ref.csv
"""
import csv
import sys
from collections import defaultdict

KT_KCAL = 0.0019872043 * 300.0          # k_B T in kcal/mol at 300 K (0.596)


def first_ok(rows, threshold, max_unvisited=0.01):
    rows = sorted(rows)                 # (t, L2, n_g, missing)
    for i, (t, *_) in enumerate(rows):
        if all(e < threshold and m <= max_unvisited * (n + m) for _, e, n, m in rows[i:]):
            return f"{t:g}"
    return "--"


def main(path):
    data = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            data[r["method"]].append((float(r["t_eff_ns"]), float(r["L2_kcalmol"]),
                                      int(r["N_g"]), int(r["missing"])))
    print(f"{'method':24s} {'kBT (0.596)':>12s} {'0.5 kcal/mol':>13s}")
    for m, rows in data.items():
        print(f"{m:24s} {first_ok(rows, KT_KCAL):>12s} {first_ok(rows, 0.5):>13s}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "l2_ref.csv")
