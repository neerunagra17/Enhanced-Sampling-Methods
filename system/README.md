# System files

- `ala_water.prmtop` — AMBER topology: Ace–Ala–Nme, ff14SB, TIP3P water, 2,767 atoms.
- `gpu_npt_final.xml` — OpenMM state after minimisation, heating to 300 K and 2 ns NPT
  equilibration at 300 K and 1 bar (positions, velocities, box vectors). This is the common
  starting state of all production runs; `common/md.py` copies positions, velocities and box
  into each NVT simulation.

Copy both files to `$PROJ/01_system/` before running (`config.sh`: `PRMTOP`, `START_XML`).
