#!/usr/bin/env python3
"""openmm_run.py - restartable OpenMM + PLUMED run (OPES, OPES+bucket) on one GPU.

Input systems: AMBER (--format amber: prmtop + rst7/inpcrd) or GROMACS (--format gromacs: topol.top + .gro).
Physics: NVT, 300 K, Langevin-middle (1/ps), PME 1.0 nm, H-bond constraints, rigid water, DT_FS (2 fs default).
Runs in chunks that end 1 step after a multiple of CHUNK steps so PLUMED's STATE file (written every CHUNK steps)
and OpenMM's state.xml are consistent; restart = re-run the same command (PLUMED gets RESTART script + STATE.rst copy).
Exit code 3 = wall-time budget used (Slurm wrapper resubmits), 0 = finished.
"""
import argparse, json, os, shutil, sys, time
import openmm as mm
from openmm import app, unit
from openmmplumed import PlumedForce

ap = argparse.ArgumentParser()
ap.add_argument("--format", choices=["amber", "gromacs"], required=True)
ap.add_argument("--top", required=True)
ap.add_argument("--coords", required=True)
ap.add_argument("--plumed-first", default="plumed_first.dat")
ap.add_argument("--plumed-restart", default="plumed_restart.dat")
ap.add_argument("--ns", type=float, default=500.0)
ap.add_argument("--chunk-ns", type=float, default=2.0)
ap.add_argument("--max-hours", type=float, default=44.0)
ap.add_argument("--dt-fs", type=float, default=2.0)
ap.add_argument("--hmr", type=float, default=0.0, help="hydrogen mass (amu) for HMR, e.g. 3.0 with dt 4 fs; 0 = off")
ap.add_argument("--seed", type=int, default=2025)
ap.add_argument("--platform", default="CUDA")
ap.add_argument("--dcd-ps", type=float, default=100.0, help="trajectory frame interval (ps); 0 = no trajectory")
a = ap.parse_args()

TEMP, FRIC = 300.0, 1.0
steps_total = int(round(a.ns * 1e6 / a.dt_fs))
steps_chunk = int(round(a.chunk_ns * 1e6 / a.dt_fs))
t_start = time.time()

kw = dict(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds,
          rigidWater=True, ewaldErrorTolerance=5e-4)
if a.hmr > 0:
    kw["hydrogenMass"] = a.hmr * unit.amu
if a.format == "amber":
    top = app.AmberPrmtopFile(a.top); crd = app.AmberInpcrdFile(a.coords)
    box = crd.boxVectors; pos = crd.positions
    system = top.createSystem(**kw)
else:
    gro = app.GromacsGroFile(a.coords)
    box = gro.getPeriodicBoxVectors(); pos = gro.positions
    top = app.GromacsTopFile(a.top, periodicBoxVectors=box)
    system = top.createSystem(**kw)

restart = os.path.exists("state.xml")
if restart:
    src = "STATE" if os.path.exists("STATE") and os.path.getsize("STATE") > 0 else "STATE.prev"
    if not os.path.exists(src) or os.path.getsize(src) == 0:
        sys.exit("restart requested but no non-empty PLUMED STATE / STATE.prev found")
    shutil.copy(src, "STATE.rst")
system.addForce(PlumedForce(open(a.plumed_restart if restart else a.plumed_first).read()))
integ = mm.LangevinMiddleIntegrator(TEMP * unit.kelvin, FRIC / unit.picosecond, a.dt_fs * unit.femtoseconds)
prog = json.load(open("progress.json")) if os.path.exists("progress.json") else {"ns_target": a.ns, "dt_fs": a.dt_fs, "parts": []}
integ.setRandomNumberSeed(a.seed + len(prog["parts"]))
platform = mm.Platform.getPlatformByName(a.platform)
sim = app.Simulation(top.topology, system, integ, platform,
                     {"Precision": "mixed"} if a.platform == "CUDA" else {})
if restart:
    sim.loadState("state.xml")
    print(f"RESTART from step {sim.context.getStepCount()}", flush=True)
else:
    sim.context.setPositions(pos)
    sim.context.setPeriodicBoxVectors(*box)
    sim.context.setVelocitiesToTemperature(TEMP * unit.kelvin, a.seed)
    sim.context.setStepCount(0)
step0 = sim.context.getStepCount()
sim.reporters.append(app.StateDataReporter("md.log", min(5000, steps_chunk), step=True, time=True, potentialEnergy=True,
                                           temperature=True, volume=True, speed=True, append=restart))
if a.dcd_ps > 0:
    sim.reporters.append(app.DCDReporter("traj.dcd", int(round(a.dcd_ps * 1000 / a.dt_fs)), append=restart and os.path.exists("traj.dcd")))
part = {"part": len(prog["parts"]) + 1, "step_start": step0, "restart": restart,
        "host": os.uname().nodename, "slurm_job": os.environ.get("SLURM_JOB_ID"), "dt_fs": a.dt_fs, "hmr": a.hmr}
prog["parts"].append(part)

while sim.context.getStepCount() < steps_total + 1:
    nxt = min((sim.context.getStepCount() // steps_chunk + 1) * steps_chunk + 1, steps_total + 1)
    n = nxt - sim.context.getStepCount()
    t0 = time.time()
    sim.step(n)
    secs = time.time() - t0
    if os.path.exists("state.xml"): shutil.copy("state.xml", "state_prev.xml")
    if os.path.exists("STATE"): shutil.copy("STATE", "STATE.prev")
    sim.saveState("state.xml")
    done = sim.context.getStepCount()
    part["step_end"] = done
    part["ns_per_day"] = round(n * a.dt_fs * 1e-6 / (secs / 86400.0), 1)
    json.dump(prog, open("progress.json", "w"), indent=1)
    print(f"step {done} ({done*a.dt_fs*1e-6:.2f} ns)  {part['ns_per_day']} ns/day", flush=True)
    if done < steps_total + 1 and a.max_hours * 3600 - (time.time() - t_start) < 1.3 * secs:
        print("wall-time budget used; exiting for resubmission", flush=True)
        sys.exit(3)
open("DONE", "w").write(f"finished {sim.context.getStepCount()} steps\n")
print("DONE", flush=True)
