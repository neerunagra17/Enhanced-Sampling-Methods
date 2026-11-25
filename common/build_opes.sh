#!/bin/bash
# build_opes.sh - build PLUMED's OPES module as a runtime library (opes/OPESmetad.so)
# for a PLUMED installation compiled WITHOUT the opes module (the conda-forge
# plumed 2.9.2 build has "module opes off").  The library is loaded in the PLUMED
# input with:   LOAD FILE=<path>/OPESmetad.so
#
# Usage:  bash build_opes.sh            (inside the activated ugp_gpu environment)
# Needs:  a C++ compiler (conda install -c conda-forge cxx-compiler, or `module load gcc`)
#         and the PLUMED headers (installed with the conda plumed package).
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

PLUMED_VERSION=$(plumed info --long-version)          # e.g. 2.9.2
echo "PLUMED $PLUMED_VERSION"
if plumed config module opes 2>/dev/null | grep -q " on"; then
  echo "This PLUMED already has the opes module: no LOAD needed."
  echo 'Set OPES_LIB="" in config.sh.'
  exit 0
fi

mkdir -p common/plumed_lib && cd common/plumed_lib
SRC=OPESmetad.cpp
if [ ! -f $SRC ]; then
  # same version as the kernel; a mismatched version will not compile/link correctly
  curl -fsSL -o $SRC "https://raw.githubusercontent.com/plumed/plumed2/v${PLUMED_VERSION}/src/opes/OPESmetad.cpp" \
    || { echo "Download failed. Copy src/opes/OPESmetad.cpp of PLUMED v${PLUMED_VERSION} here by hand."; exit 1; }
fi

# PLUMED's own `plumed mklib` hard-codes the compiler of the conda build machine,
# which does not exist on the cluster. Re-use its flags with a local compiler.
ROOT=$(plumed info --root)
OPTS="$ROOT/src/config/compile_options.sh"
CXX_LOCAL=${CXX:-$(command -v x86_64-conda-linux-gnu-c++ || command -v g++ || command -v c++)}
[ -n "$CXX_LOCAL" ] || { echo "No C++ compiler found"; exit 1; }
echo "compiler: $CXX_LOCAL"
compile=$(bash -c "source $OPTS; echo \"\$compile\"")
link=$(bash -c "source $OPTS; echo \"\$link_installed\"")
compile="$CXX_LOCAL ${compile#* }"
link="$CXX_LOCAL ${link#* }"
export PLUMED_INCLUDEDIR=$(plumed info --include-dir)
export PLUMED_PROGRAM_NAME=plumed
rm -f OPESmetad.o OPESmetad.so
eval "$compile OPESmetad.o $SRC"
# kernel library listed AFTER the object (with --no-as-needed) so that OPESmetad.so records
# its dependency on libplumedKernel.so; otherwise LOAD fails inside openmm-plumed with
# "undefined symbol: _ZTIN4PLMD6ActionE" (it works in `plumed driver` only by accident).
KERNEL_LIB=${PLUMED_KERNEL:-$(cd "$ROOT/.." && pwd)/libplumedKernel.so}
[ -f "$KERNEL_LIB" ] || { echo "cannot find libplumedKernel.so (set PLUMED_KERNEL)"; exit 1; }
eval "$link OPESmetad.so OPESmetad.o -Wl,--no-as-needed $KERNEL_LIB"
ldd OPESmetad.so | grep -q plumedKernel || { echo "OPESmetad.so is not linked to the PLUMED kernel"; exit 1; }
echo "built $(pwd)/OPESmetad.so"

# test: the action must be recognised after LOAD
cat > test_opes.dat << EOT
LOAD FILE=$(pwd)/OPESmetad.so
d: DISTANCE ATOMS=1,2
o: OPES_METAD ARG=d PACE=10 BARRIER=10 TEMP=300
EOT
printf "2\n1 1 1\nX 0 0 0\nX 0 0 1\n" > test.xyz
rm -f KERNELS
plumed driver --plumed test_opes.dat --ixyz test.xyz --kt 2.494 > test_opes.log 2>&1 \
  && grep -q "Action OPES_METAD" test_opes.log && echo "OK: OPES_METAD loads and runs" \
  || { echo "FAILED: see common/plumed_lib/test_opes.log"; exit 1; }
rm -f KERNELS test.xyz
echo "OPES inputs will use:  LOAD FILE=$(pwd)/OPESmetad.so"
