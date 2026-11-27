import sys, openmm as mm
from openmm import unit
plat = sys.argv[1]
s = mm.System(); nb = mm.NonbondedForce(); nb.setNonbondedMethod(nb.PME)
s.setDefaultPeriodicBoxVectors(mm.Vec3(3,0,0), mm.Vec3(0,3,0), mm.Vec3(0,0,3))
import random
for i in range(200):
    s.addParticle(12.0); nb.addParticle(0.1 if i%2 else -0.1, 0.3, 0.5)
s.addForce(nb)
ctx = mm.Context(s, mm.LangevinMiddleIntegrator(300, 1, 0.002), mm.Platform.getPlatformByName(plat), {"Precision": "mixed"})
ctx.setPositions([mm.Vec3(random.random()*3, random.random()*3, random.random()*3) for _ in range(200)])
ctx.getIntegrator().step(10); ctx.getState(getEnergy=True)
print("context OK on", plat, flush=True)
