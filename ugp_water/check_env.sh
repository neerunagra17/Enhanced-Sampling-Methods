#!/bin/bash
# check_env.sh - run once on the cluster before starting.
# Prints what is available; fix anything marked MISSING before running.

echo "== AMBER =="
[ -n "$AMBERHOME" ] && echo "AMBERHOME = $AMBERHOME" || echo "MISSING: AMBERHOME not set (source amber.sh)"
for exe in tleap sander sander.MPI pmemd pmemd.cuda cpptraj; do
  command -v $exe >/dev/null && echo "  found   $exe" || echo "  missing $exe"
done

echo "== PLUMED =="
if command -v plumed >/dev/null; then
  plumed info --version
  # OPES needs the 'opes' module; the lab's TAMD/TASS needs EXTENDED_LAGRANGIAN
  for action in METAD OPES_METAD RESTRAINT UPPER_WALLS EXTENDED_LAGRANGIAN; do
    plumed manual --action $action >/dev/null 2>&1 \
      && echo "  found   $action" || echo "  MISSING $action"
  done
else
  echo "MISSING: plumed executable (module load plumed?)"
fi
[ -n "$PLUMED_KERNEL" ] && echo "PLUMED_KERNEL = $PLUMED_KERNEL" \
  || echo "note: PLUMED_KERNEL not set (sander needs it to load PLUMED at run time)"

echo "== Python =="
python3 -c "import numpy, scipy, matplotlib; print('numpy', numpy.__version__, '| scipy', scipy.__version__, '| matplotlib', matplotlib.__version__)" \
  || echo "MISSING: pip install --user numpy scipy matplotlib"
