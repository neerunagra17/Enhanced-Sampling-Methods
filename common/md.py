#!/usr/bin/env python3
"""md.py - NVT production/preparation run with OpenMM + PLUMED (openmm-plumed),
with checkpoint/restart. Used by every method in this project.

    python md.py --plumed plumed.dat --outdir RUN --ns 15.625 --start START.xml

Behaviour
  * The system is rebuilt from the AMBER prmtop with the SAME settings as the
    equilibration (PME, HBonds constraints, rigid water, 2 fs, Langevin 300 K)
    but WITHOUT a barostat (NVT).
  * Fresh start: positions, velocities and box vectors are copied one by one
    from --start (an OpenMM State XML, e.g. the NPT result gpu_npt_final.xml).
    Context.setState()/Simulation.loadState() is NOT used, because the NPT
    state carries the barostat parameters (MonteCarloPressure, ...) that do
    not exist in the NVT context. Time and step count are reset to 0.
  * Restart: if RUN/state.xml exists, the run continues from it (positions,
    velocities, box, time, step count), PLUMED is started in RESTART mode, and
    PLUMED output files are first cut back to the checkpoint time so that no
    bias or COLVAR data from beyond the checkpoint is kept (see sync_plumed_files).
  * Checkpoints (RUN/state.xml, portable XML, written atomically) every --chk-ps.
    Binary .chk files are not used: they are only valid on the same GPU/driver/
    OpenMM build, while state.xml can be restarted on any node.
  * The run stops cleanly at a checkpoint when the next chunk would exceed
    --max-hours. Exit code 0 = finished (RUN/DONE written), 3 = stopped for
    wall time (resubmit), anything else = error.
"""
import argparse
import json
import os
import re
import shutil
import sys
import time

import openmm as mm
import openmm.app as app
import openmm.unit as u
from openmmplumed import PlumedForce

PLUMED_GROUP = 30          # force group of the PLUMED bias
EXIT_DONE, EXIT_WALLTIME = 0, 3
# PLUMED output files whose first column is time (ps) and that are cut back on restart
SYNC_PREFIXES = ("COLVAR", "HILLS", "KERNELS")
# PLUMED files that are complete snapshots (OPES STATE); copied to <name>.chk at every
# checkpoint and restored from there on restart
SNAP_PREFIXES = ("STATE",)


def log(msg):
    print(time.strftime("[%Y-%m-%d %H:%M:%S] ") + msg, flush=True)


def parse():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    e = os.environ.get
    ap.add_argument("--prmtop", default=e("PRMTOP"))
    ap.add_argument("--start", help="State XML for a fresh start (ignored when restarting)")
    ap.add_argument("--plumed", required=True, help="PLUMED input file")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ns", type=float, required=True, help="total length of this run (ns)")
    ap.add_argument("--max-hours", type=float, default=1e9,
                    help="stop cleanly before this many hours (wall time minus margin)")
    ap.add_argument("--platform", default=e("PLATFORM", "OpenCL"))
    ap.add_argument("--precision", default=e("PRECISION", "mixed"))
    ap.add_argument("--temp", type=float, default=float(e("TEMP", 300)))
    ap.add_argument("--dt-fs", type=float, default=float(e("DT_FS", 2)))
    ap.add_argument("--friction", type=float, default=float(e("FRICTION", 1.0)))
    ap.add_argument("--cutoff-nm", type=float, default=float(e("CUTOFF_NM", 1.0)))
    ap.add_argument("--ewald-tol", type=float, default=float(e("EWALD_TOL", 5e-4)))
    ap.add_argument("--chk-ps", type=float, default=float(e("CHK_PS", 100)))
    ap.add_argument("--report-ps", type=float, default=float(e("REPORT_PS", 10)))
    ap.add_argument("--dcd-ps", type=float, default=float(e("DCD_PS", 10)))
    ap.add_argument("--seed", type=int, default=0, help="Langevin seed (0 = random)")
    ap.add_argument("--stop-after-chunks", type=int, default=0,
                    help="testing only: stop (exit 3) after this many checkpoint chunks")
    return ap.parse_args()


# ----------------------------------------------------------------------------
def build_system(a):
    prmtop = app.AmberPrmtopFile(a.prmtop)
    system = prmtop.createSystem(nonbondedMethod=app.PME,
                                 nonbondedCutoff=a.cutoff_nm * u.nanometer,
                                 constraints=app.HBonds, rigidWater=True,
                                 ewaldErrorTolerance=a.ewald_tol)
    for f in system.getForces():
        if isinstance(f, (mm.MonteCarloBarostat, mm.MonteCarloAnisotropicBarostat,
                          mm.MonteCarloFlexibleBarostat)):
            raise RuntimeError("barostat found in the NVT system")
    return prmtop, system


def make_platform(a):
    plat = mm.Platform.getPlatformByName(a.platform)
    props = {}
    if a.platform in ("CUDA", "OpenCL"):
        props["Precision"] = a.precision
    return plat, props


def read_state(path):
    with open(path) as fh:
        return mm.XmlSerializer.deserialize(fh.read())


def transfer_state(context, st, natoms, keep_clock):
    """Copy box, positions and velocities from a State (NPT or NVT) into the context."""
    pos = st.getPositions(asNumpy=True)
    if len(pos) != natoms:
        raise RuntimeError(f"state has {len(pos)} atoms, topology has {natoms}")
    try:
        vel = st.getVelocities(asNumpy=True)
    except mm.OpenMMException:
        raise RuntimeError("state XML has no velocities; save it with "
                           "getState(getPositions=True, getVelocities=True)")
    context.setPeriodicBoxVectors(*st.getPeriodicBoxVectors())
    context.setPositions(pos)
    context.setVelocities(vel)
    context.applyConstraints(1e-6)
    context.applyVelocityConstraints(1e-6)
    if keep_clock:
        context.setTime(st.getTime())
        context.setStepCount(st.getStepCount())
    else:
        context.setTime(0.0)
        context.setStepCount(0)


def save_state(context, outdir):
    """Atomic write of the checkpoint (state.xml; previous one kept as state_prev.xml).

    The energy call first makes PLUMED process the current step (calc + update:
    COLVAR row, hill/kernel deposition, OPES STATE written when the step is a
    multiple of STATE_WSTRIDE). openmm-plumed runs update() only once per step,
    so the next MD step does not repeat it. PLUMED snapshot files are then copied
    so that the restart files always match state.xml.
    """
    context.getState(getEnergy=True, groups={PLUMED_GROUP})
    for name in os.listdir(outdir):
        if name.startswith(SNAP_PREFIXES) and not name.endswith((".chk", ".tmp", ".bak")):
            src = os.path.join(outdir, name)
            shutil.copyfile(src, src + ".chk.tmp")
            os.replace(src + ".chk.tmp", src + ".chk")
    st = context.getState(getPositions=True, getVelocities=True)
    tmp = os.path.join(outdir, "state.xml.tmp")
    with open(tmp, "w") as fh:
        fh.write(mm.XmlSerializer.serialize(st))
    cur = os.path.join(outdir, "state.xml")
    if os.path.exists(cur):
        shutil.copyfile(cur, os.path.join(outdir, "state_prev.xml"))
    os.replace(tmp, cur)
    return st


def sync_plumed_files(outdir, t_chk_ps, dt_ps):
    """Cut PLUMED time-series files back to the checkpoint time.

    After a hard kill (node failure, scancel), PLUMED may have written data
    beyond the last checkpoint; keeping it would double-count hills/kernels.
    Rows with time > t_chk are removed. The row at exactly t_chk is removed
    from COLVAR files (it is printed again when the run restarts at that step)
    and kept in HILLS/KERNELS (PLUMED does not deposit on the first step of a
    restarted run). Header lines are kept.
    """
    eps = 0.25 * dt_ps
    for name in sorted(os.listdir(outdir)):
        if not name.startswith(SYNC_PREFIXES) or name.endswith((".bak", ".tmp")):
            continue
        path = os.path.join(outdir, name)
        limit = t_chk_ps - eps if name.startswith("COLVAR") else t_chk_ps + eps
        kept, dropped = [], 0
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    kept.append(line)
                    continue
                try:
                    t = float(line.split()[0])
                except ValueError:
                    kept.append(line)
                    continue
                if t <= limit:
                    kept.append(line)
                else:
                    dropped += 1
        if dropped:
            shutil.copyfile(path, path + ".bak")
            with open(path + ".tmp", "w") as fh:
                fh.writelines(kept)
            os.replace(path + ".tmp", path)
            log(f"restart sync: removed {dropped} rows after t={t_chk_ps:.1f} ps from {name}")


# ----------------------------------------------------------------------------
def main():
    a = parse()
    t_wall0 = time.time()
    outdir = os.path.abspath(a.outdir)
    os.makedirs(outdir, exist_ok=True)
    plumed_in = os.path.abspath(a.plumed)
    prmtop_path = os.path.abspath(a.prmtop)
    start_path = os.path.abspath(a.start) if a.start else None
    a.prmtop = prmtop_path

    done_flag = os.path.join(outdir, "DONE")
    if os.path.exists(done_flag):
        log(f"{outdir}: already finished (DONE exists)")
        return EXIT_DONE

    dt_ps = a.dt_fs / 1000.0
    nsteps_total = int(round(a.ns * 1000.0 / dt_ps))
    chk_steps = max(1, int(round(a.chk_ps / dt_ps)))
    rep_steps = max(1, int(round(a.report_ps / dt_ps)))
    dcd_steps = int(round(a.dcd_ps / dt_ps)) if a.dcd_ps > 0 else 0
    if chk_steps % rep_steps:
        raise SystemExit("--chk-ps must be a multiple of --report-ps")

    state_xml = os.path.join(outdir, "state.xml")
    restart = os.path.exists(state_xml)
    prog_path = os.path.join(outdir, "progress.json")
    prog = json.load(open(prog_path)) if os.path.exists(prog_path) else {
        "ns_target": a.ns, "dt_fs": a.dt_fs, "start": start_path, "parts": []}

    prmtop, system = build_system(a)
    natoms = system.getNumParticles()

    os.chdir(outdir)                      # PLUMED reads/writes relative to here
    if restart:
        st0 = read_state(state_xml)
        sync_plumed_files(outdir, st0.getTime().value_in_unit(u.picosecond), dt_ps)
        for name in os.listdir(outdir):           # restore OPES STATE matching state.xml
            if name.startswith(SNAP_PREFIXES) and name.endswith(".chk"):
                shutil.copyfile(name, name[:-4])
                log(f"restart: restored {name[:-4]} from {name}")
    else:
        if not start_path:
            raise SystemExit("fresh start needs --start")
        st0 = read_state(start_path)
        shutil.copyfile(plumed_in, os.path.join(outdir, "plumed_input_used.dat"))

    with open(plumed_in) as fh:
        script = fh.read()
    # OPES STATE must be written exactly at md.py checkpoints (see save_state)
    for m_ in re.finditer(r"STATE_WSTRIDE=(\d+)", script):
        if int(m_.group(1)) != chk_steps:
            raise SystemExit(f"STATE_WSTRIDE={m_.group(1)} in {plumed_in} differs from the checkpoint "
                             f"interval ({chk_steps} steps = --chk-ps {a.chk_ps}); regenerate plumed.dat")
    pf = PlumedForce(script)
    pf.setTemperature(a.temp * u.kelvin)   # kBT for actions without TEMP=
    pf.setRestart(restart)                 # PLUMED RESTART mode (append, re-read HILLS/STATE)
    pf.setForceGroup(PLUMED_GROUP)
    system.addForce(pf)

    integ = mm.LangevinMiddleIntegrator(a.temp * u.kelvin, a.friction / u.picosecond,
                                        a.dt_fs * u.femtosecond)
    integ.setRandomNumberSeed(a.seed)
    plat, props = make_platform(a)
    sim = app.Simulation(prmtop.topology, system, integ, plat, props)
    transfer_state(sim.context, st0, natoms, keep_clock=restart)
    step0 = sim.context.getStepCount()
    log(f"{'RESTART' if restart else 'START'} {outdir}: step {step0}/{nsteps_total} "
        f"({step0*dt_ps/1000:.3f}/{a.ns} ns), platform {plat.getName()} {props}, "
        f"{natoms} atoms, PLUMED restart={restart}")
    if step0 >= nsteps_total:
        open(done_flag, "w").write("finished\n")
        return EXIT_DONE
    if a.max_hours < 0.05:
        log(f"only {a.max_hours*60:.1f} min of wall time left: not starting")
        return EXIT_WALLTIME

    part = len(prog["parts"]) + 1
    sim.reporters.append(app.StateDataReporter(
        os.path.join(outdir, "md.log"), rep_steps, step=True, time=True,
        potentialEnergy=True, temperature=True, volume=True, speed=True,
        append=restart and os.path.exists(os.path.join(outdir, "md.log"))))
    if dcd_steps:
        # one trajectory file per part; frames of a killed part after its last
        # checkpoint are listed as invalid in progress.json
        sim.reporters.append(app.DCDReporter(os.path.join(outdir, f"traj_{part:03d}.dcd"),
                                             dcd_steps))
    prog["parts"].append({"part": part, "step_start": step0, "step_end": step0,
                          "restart": restart, "host": os.uname().nodename,
                          "slurm_job": os.environ.get("SLURM_JOB_ID", "")})

    t_chunk = None
    nchunks = 0
    step = step0
    status = EXIT_DONE
    while step < nsteps_total:
        n = min(chk_steps - step % chk_steps, nsteps_total - step)
        c0 = time.time()
        sim.step(n)
        step = sim.context.getStepCount()
        save_state(sim.context, outdir)
        dt_wall = time.time() - c0
        t_chunk = dt_wall * chk_steps / n
        nsday = n * dt_ps / 1000.0 / (dt_wall / 86400.0)
        prog["parts"][-1].update(step_end=step, ns_per_day=round(nsday, 1))
        json.dump(prog, open(prog_path + ".tmp", "w"), indent=1)
        os.replace(prog_path + ".tmp", prog_path)
        elapsed_h = (time.time() - t_wall0) / 3600.0
        nchunks += 1
        if step < nsteps_total and a.stop_after_chunks and nchunks >= a.stop_after_chunks:
            log(f"test stop after {nchunks} chunks at step {step}")
            status = EXIT_WALLTIME
            break
        if step < nsteps_total and elapsed_h + 1.2 * t_chunk / 3600.0 > a.max_hours:
            log(f"wall-time limit: stopping at step {step} ({step*dt_ps/1000:.3f} ns), "
                f"{nsday:.0f} ns/day")
            status = EXIT_WALLTIME
            break
    if step >= nsteps_total:
        open(done_flag, "w").write(f"finished {step} steps\n")
        log(f"finished {outdir}: {step} steps = {step*dt_ps/1000:.3f} ns")
    del sim                                # finalises PLUMED (flushes its files)
    return status


if __name__ == "__main__":
    sys.exit(main())
