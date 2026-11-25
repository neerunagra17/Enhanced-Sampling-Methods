#!/usr/bin/env python3
"""check_outputs.py DIR - check the smoke/restart test written by slurm/job_check.sh.

For each method folder in DIR (us, metad, opes, tamd, tass):
  * the run finished in two parts, the second one a restart
  * COLVAR time is strictly increasing with a constant stride (no gap, no
    duplicate at the restart) and reaches the requested length
  * MetaD/TASS: HILLS strictly increasing, one hill per PACE, none lost or
    doubled at the restart, and PLUMED reported re-reading HILLS
  * OPES: PLUMED restarted from STATE; KERNELS strictly increasing
  * TAMD/TASS: size of the auxiliary-variable reset at the restart
  * speed (ns/day) and the GPU time needed for 500 ns
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
import fesutils as fu  # noqa: E402

nfail = 0


def say(level, msg):
    global nfail
    nfail += level == "FAIL"
    print(f"[{level:4s}] {msg}")


def series(path):
    """Raw first column (time) of a PLUMED file, without any filtering."""
    t = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            t.append(float(line.split()[0]))
    return np.array(t)


def check_times(label, t, stride, t_end=None, t_first=None):
    d = np.diff(t)
    ok = len(t) > 2 and np.all(d > 0) and np.allclose(d, stride, atol=1e-6)
    msg = f"{label}: {len(t)} rows, {t[0]:g}..{t[-1]:g} ps, stride {np.median(d):g} ps"
    if not ok:
        bad = np.where(~np.isclose(d, stride, atol=1e-6))[0][:5]
        msg += f"; irregular at rows {bad.tolist()} (times {[t[i] for i in bad]}->{[t[i+1] for i in bad]})"
    if t_end is not None and abs(t[-1] - t_end) > 1e-6:
        ok = False
        msg += f"; expected last time {t_end:g}"
    if t_first is not None and abs(t[0] - t_first) > 1e-6:
        ok = False
        msg += f"; expected first time {t_first:g}"
    say("OK" if ok else "FAIL", msg)


def main(top):
    speeds = {}
    for m in ("us", "metad", "opes", "tamd", "tass"):
        d = os.path.join(top, m)
        print(f"\n== {m}")
        if not os.path.isdir(d):
            say("FAIL", "folder missing")
            continue
        prog = json.load(open(os.path.join(d, "progress.json")))
        parts = prog["parts"]
        ns_total = prog["ns_target"] * 1000
        ok = os.path.exists(os.path.join(d, "DONE")) and len(parts) == 2 and parts[1]["restart"]
        say("OK" if ok else "FAIL", f"finished in {len(parts)} parts, second is restart: "
            f"{len(parts) == 2 and parts[-1]['restart']}")
        t_restart = parts[0]["step_end"] * prog["dt_fs"] / 1000.0
        speeds[m] = parts[-1].get("ns_per_day", np.nan)
        log2 = open(os.path.join(d, "part2.log")).read() if os.path.exists(os.path.join(d, "part2.log")) else ""

        stride = 0.1 if m in ("us", "tass") else 1.0
        t = series(os.path.join(d, "COLVAR"))
        check_times("COLVAR", t, stride, t_end=ns_total, t_first=0.0)
        cv = fu.load_colvar(os.path.join(d, "COLVAR"))

        if m in ("metad", "tass"):
            h = series(os.path.join(d, "HILLS"))
            check_times("HILLS", h, 1.0, t_end=ns_total, t_first=1.0)
            say("OK" if "Restarting from HILLS" in log2 else "FAIL",
                "PLUMED re-read HILLS on restart")
        if m == "opes":
            say("OK" if "restarting from: STATE" in log2 else "FAIL", "PLUMED restarted OPES from STATE")
            k = series(os.path.join(d, "KERNELS"))
            if len(k) > 2:
                check_times("KERNELS", k, 1.0, t_end=ns_total)
            else:
                say("WARN", f"only {len(k)} kernels (adaptive sigma skips the first 10)")
            say("OK" if os.path.exists(os.path.join(d, "STATE.chk")) else "FAIL", "STATE.chk snapshot present")
        if m in ("tamd", "tass"):
            diff = fu.wrap(cv["ext.phi_fict"] - cv["phi"])
            i = int(np.argmin(np.abs(cv["time"] - t_restart)))
            rms = np.sqrt(np.mean(diff ** 2))
            say("INFO", f"s - q (phi): rms {rms:.3f} rad; at the restart (t={cv['time'][i]:g} ps) "
                        f"{diff[i]:+.3f} rad (EXTENDED_LAGRANGIAN resets s=q, v_s=0 on restart)")
        if m == "us":
            say("INFO", f"phi mean {np.mean(cv['phi']):.3f} rad (restraint at -1.4)")

    print("\n== speed (last part) and GPU time for production")
    for m, v in speeds.items():
        if np.isfinite(v) and v > 0:
            print(f"  {m:6s} {v:8.1f} ns/day  ->  500 ns needs {500 / v:6.1f} GPU-days"
                  + ("  (32 windows run in parallel)" if m in ("us", "tass") else ""))
    print("  (short test runs with 5 ps checkpoints: production speed is usually a little higher)")
    print(f"\n{nfail} FAIL")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
