# windows.sh - window centres (rad), in the ORDER used for the chained preparation.
# 32 centres, 0.2 rad apart. Start near the equilibrated structure (phi ~ -1.4 in
# water is typical; check the first frame) and walk outwards in both directions.
# Because phi is periodic, -3.2 rad is the same as +3.083 rad, so the second
# branch starts from the restart of the -3.2 window.
BRANCH1="-1.4 -1.6 -1.8 -2.0 -2.2 -2.4 -2.6 -2.8 -3.0 -3.2"
BRANCH2="3.0 2.8 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 0.0 -0.2 -0.4 -0.6 -0.8 -1.0 -1.2"
ALL="$BRANCH1 $BRANCH2"
export KAPPA=500.0
