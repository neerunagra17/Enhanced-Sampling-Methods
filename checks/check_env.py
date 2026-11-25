#!/usr/bin/env python3
"""check_env.py - verify the software and the system before any production run.

Run on a GPU node (submit.sh check does this). Every line is OK / WARN / FAIL.
  1. versions: OpenMM, openmm-plumed, PLUMED kernel; PLUMED modules and actions used here
  2. platforms: OpenCL present, which device is used, CUDA_VISIBLE_DEVICES
  3. topology: phi/psi atom indices 5,7,9,15 / 7,9,15,17 are C(ACE) N CA C(ALA) N(NME)
  4. start state: atoms, box, velocities present, barostat parameters (why loadState fails)
  5. PlumedForce on the GPU: PLUMED forces and energy on OpenCL agree with the CPU platform
"""
import os
import subprocess
import sys
from importlib import metadata

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u

status = {"FAIL": 0, "WARN": 0}


def say(level, msg):
    if level in status:
        status[level] += 1
    print(f"[{level:4s}] {msg}", flush=True)


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.returncode, (r.stdout + r.stderr).strip()


env = os.environ
prmtop_path, start_path = env["PRMTOP"], env["START_XML"]
platform_name, opes_lib = env.get("PLATFORM", "OpenCL"), env.get("OPES_LIB", "")

# 1 ---------------------------------------------------------------- versions
say("INFO", f"OpenMM {mm.__version__}")
for pkg in ("openmm-plumed", "openmmplumed"):
    try:
        say("INFO", f"{pkg} {metadata.version(pkg)}")
        break
    except metadata.PackageNotFoundError:
        pass
else:
    rc, out = sh("conda list openmm-plumed | tail -1")
    say("INFO", f"openmm-plumed (conda): {out}")
try:
    from openmmplumed import PlumedForce
    say("OK", "openmmplumed imports")
except Exception as e:
    say("FAIL", f"cannot import openmmplumed: {e}")
    sys.exit(1)
rc, ver = sh("plumed info --long-version")
say("OK" if rc == 0 else "FAIL", f"PLUMED executable: {ver}")
rc, out = sh("plumed config show | grep '^module'")
mods = {l.split()[1]: l.split()[2] for l in out.splitlines() if len(l.split()) > 2}
for m in ("bias", "colvar", "function", "opes"):
    say("INFO", f"PLUMED module {m}: {mods.get(m, '?')}")
for act in ("TORSION", "RESTRAINT", "METAD", "EXTENDED_LAGRANGIAN"):
    rc, _ = sh(f"plumed manual --action {act}")
    say("OK" if rc == 0 else "FAIL", f"PLUMED action {act}")
if mods.get("opes") == "on":
    say("OK", "OPES module built in (set OPES_LIB=\"\" in config.sh)")
elif opes_lib and os.path.isfile(opes_lib):
    say("OK", f"OPES via LOAD: {opes_lib}")
else:
    say("FAIL", "OPES_METAD unavailable: run  bash common/build_opes.sh")

# 2 ---------------------------------------------------------------- platforms
names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
say("INFO", f"OpenMM platforms: {names}")
say("OK" if platform_name in names else "FAIL", f"requested platform {platform_name}")
say("INFO", f"CUDA_VISIBLE_DEVICES={env.get('CUDA_VISIBLE_DEVICES', 'unset')} "
            f"SLURM_JOB_GPUS={env.get('SLURM_JOB_GPUS', 'unset')}")

# 3 ---------------------------------------------------------------- topology
prmtop = app.AmberPrmtopFile(prmtop_path)
atoms = list(prmtop.topology.atoms())
expect = {5: ("C", "ACE"), 7: ("N", "ALA"), 9: ("CA", "ALA"), 15: ("C", "ALA"), 17: ("N", "NME")}
for i, (name, res) in expect.items():
    a = atoms[i - 1]
    ok = a.name == name and a.residue.name == res
    say("OK" if ok else "FAIL", f"PLUMED atom {i}: {a.name} {a.residue.name} (expected {name} {res})")
say("INFO", f"{len(atoms)} atoms, {prmtop.topology.getNumResidues()} residues")

# 4 ---------------------------------------------------------------- start state
st = mm.XmlSerializer.deserialize(open(start_path).read())
pos = st.getPositions(asNumpy=True)
say("OK" if len(pos) == len(atoms) else "FAIL", f"start state atoms {len(pos)}")
try:
    vel = st.getVelocities(asNumpy=True).value_in_unit(u.nanometer / u.picosecond)
    say("OK", f"start state has velocities (rms {np.sqrt((vel**2).sum(1).mean()):.3f} nm/ps)")
except mm.OpenMMException:
    say("FAIL", "start state has no velocities")
box = st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer)
vol = abs(np.linalg.det(box))
say("INFO", f"box vectors (nm): {np.round(box, 4).tolist()}  volume {vol:.3f} nm^3")
say("INFO", f"start state time {st.getTime()}, step {st.getStepCount()}")
try:
    params = st.getParameters()
    say("INFO", f"start state parameters {dict(params)}: Simulation.loadState() would fail "
                "in the NVT system; md.py copies positions/velocities/box instead")
except mm.OpenMMException:
    say("INFO", "start state has no parameters")

# 5 ---------------------------------------------------------------- PlumedForce on GPU vs CPU
def plumed_forces(platform, props):
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=float(env.get("CUTOFF_NM", 1.0)) * u.nanometer,
                                 constraints=app.HBonds, rigidWater=True)
    script = ("phi: TORSION ATOMS=5,7,9,15\npsi: TORSION ATOMS=7,9,15,17\n"
              "RESTRAINT ARG=phi,psi AT=0.5,0.5 KAPPA=500,500\n")
    pf = PlumedForce(script)
    pf.setForceGroup(30)
    pf.setTemperature(300 * u.kelvin)
    system.addForce(pf)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName(platform), props)
    ctx.setPeriodicBoxVectors(*st.getPeriodicBoxVectors())
    ctx.setPositions(pos)
    s = ctx.getState(getForces=True, getEnergy=True, groups={30})
    f = s.getForces(asNumpy=True).value_in_unit(u.kilojoule_per_mole / u.nanometer)
    e = s.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    dev = ctx.getPlatform().getPropertyValue(ctx, "DeviceName") if platform in ("OpenCL", "CUDA") else "cpu"
    del ctx
    return e, f, dev


try:
    e_ref, f_ref, _ = plumed_forces("CPU" if "CPU" in names else "Reference", {})
    if platform_name in ("OpenCL", "CUDA"):
        e_gpu, f_gpu, dev = plumed_forces(platform_name, {"Precision": env.get("PRECISION", "mixed")})
        say("INFO", f"{platform_name} device: {dev}")
        df = np.abs(f_gpu - f_ref).max()
        ok = abs(e_gpu - e_ref) < 1e-3 * max(1.0, abs(e_ref)) and df < 1e-2 * max(1.0, np.abs(f_ref).max())
        say("OK" if ok else "FAIL", f"PLUMED bias on {platform_name}: E={e_gpu:.4f} vs CPU {e_ref:.4f} kJ/mol, "
                                     f"max |dF|={df:.2e} kJ/mol/nm")
    else:
        say("WARN", f"platform {platform_name}: GPU comparison skipped (E_PLUMED={e_ref:.4f})")
except Exception as e:
    say("FAIL", f"PlumedForce test failed: {e}")

print(f"\nsummary: {status['FAIL']} FAIL, {status['WARN']} WARN")
sys.exit(1 if status["FAIL"] else 0)
