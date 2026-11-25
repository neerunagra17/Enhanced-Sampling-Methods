#!/usr/bin/env python3
"""status.py - progress of every run (reads progress.json written by md.py)."""
import glob
import json
import os

E = os.environ


def row(d):
    p = os.path.join(d, "progress.json")
    if not os.path.exists(p):
        return None
    g = json.load(open(p))
    parts = g["parts"]
    done_ns = parts[-1]["step_end"] * g["dt_fs"] / 1e6 if parts else 0.0
    speed = parts[-1].get("ns_per_day", 0) if parts else 0
    return done_ns, g["ns_target"], speed, len(parts), os.path.exists(os.path.join(d, "DONE"))


def show(title, dirs):
    rows = [(d, row(d)) for d in dirs]
    rows = [(d, r) for d, r in rows if r]
    if not rows:
        print(f"{title:22s} not started")
        return
    tot = sum(r[0] for _, r in rows)
    tgt = sum(r[1] for _, r in rows)
    ndone = sum(r[4] for _, r in rows)
    print(f"{title:22s} {tot:9.3f} / {tgt:7.3f} ns   runs finished {ndone}/{len(dirs)}")
    for d, (ns, t, sp, nparts, done) in rows:
        if not done or len(rows) == 1:
            print(f"   {os.path.basename(d):14s} {ns:9.3f}/{t:g} ns  {sp:7.1f} ns/day  parts {nparts}"
                  f"{'  DONE' if done else ''}")


win = lambda base: sorted(glob.glob(os.path.join(base, "win_*")), key=lambda p: float(p.split("win_")[1]))
for name, key in (("US prep", "US_DIR"), ("TASS prep", "TASS_DIR")):
    show(name, win(os.path.join(E[key], "prep")))
show("US production", win(E["US_DIR"]))
show("TASS production", win(E["TASS_DIR"]))
for name, key in (("WT-MetaD", "METAD_DIR"), ("TAMD", "TAMD_DIR")):
    show(name, [E[key]])
