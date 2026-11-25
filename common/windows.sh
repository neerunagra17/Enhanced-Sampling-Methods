# windows.sh - umbrella/TASS window centres (rad).
# Order = order of the chained 100 ps preparation. 32 centres, 0.2 rad apart, starting
# near the equilibrated structure and walking outwards in both directions. phi is periodic:
# -3.2 rad is the same point as +3.083 rad, so branch 2 starts from the -3.2 window.
# (Because 32 x 0.2 = 6.4 > 2*pi, the windows at -3.2 and 3.0 are only 0.083 rad apart;
#  WHAM handles unequal spacing, this only adds sampling near +-pi.)
BRANCH1="-1.4 -1.6 -1.8 -2.0 -2.2 -2.4 -2.6 -2.8 -3.0 -3.2"
BRANCH2="3.0 2.8 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 0.0 -0.2 -0.4 -0.6 -0.8 -1.0 -1.2"
BRANCH2_FROM="-3.2"
ALL="$BRANCH1 $BRANCH2"
KAPPA_US=500.0          # kJ/mol/rad^2, umbrella on phi (US) and on aux. phi (TASS)
# window centre for array index i (0..31)
centre_of () { echo $ALL | awk -v i=$(( $1 + 1 )) '{print $i}'; }
