#!/usr/bin/env python3
"""metad_checkpoints.py - MetaD FES at several t_eff with PLUMED sum_hills.

For every t in --teff (ns) the hills deposited before t are written to a
temporary file and summed with 'plumed sum_hills', giving
    fes_metad_<t>ns.dat   (phi psi F ..., kJ/mol; grid from sum_hills)
read_fes() in fesutils puts these on the common grid automatically.

    python3 metad_checkpoints.py --hills HILLS --teff 10,20,50,100 --prefix fes_metad

For well-tempered MetaD, sum_hills already returns the free energy
(F = -(gamma/(gamma-1)) V, Eq. WT-MetaD of the report); no extra scaling.
"""
import argparse
import os
import subprocess
import tempfile
import numpy as np
import fesutils as fu


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hills", default="HILLS")
    ap.add_argument("--teff", required=True, help="comma list, ns")
    ap.add_argument("--prefix", default="fes_metad")
    ap.add_argument("--nbins", type=int, default=fu.NBINS)
    ap.add_argument("--plumed", default="plumed", help="plumed executable")
    a = ap.parse_args()

    header, lines = [], []
    with open(a.hills) as fh:
        for line in fh:
            if line.startswith("#"):
                if not lines:                     # keep only the first header block
                    header.append(line)
            elif line.strip():
                lines.append(line)
    t = fu.continuous_time([float(l.split()[0]) for l in lines])
    print(f"{a.hills}: {len(lines)} hills, last at {t[-1]/1000:.2f} ns")

    nb = a.nbins
    for t_ns in fu.parse_list(a.teff):
        n = int(np.searchsorted(t, t_ns * 1000.0, side="right"))
        if t_ns * 1000.0 > t[-1] + (t[1] - t[0]):
            print(f"  skip {t_ns:g} ns (run too short)")
            continue
        with tempfile.NamedTemporaryFile("w", suffix=".hills", delete=False) as tmp:
            tmp.writelines(header)
            tmp.writelines(lines[:n])
        out = f"{a.prefix}_{fu.fmt_ns(t_ns)}ns.dat"
        cmd = [a.plumed, "sum_hills", "--hills", tmp.name, "--outfile", out,
               "--mintozero", "--min", "-pi,-pi", "--max", "pi,pi",
               "--bin", f"{nb},{nb}"]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        os.unlink(tmp.name)
        print(f"  wrote {out}  ({n} hills)")


if __name__ == "__main__":
    main()
