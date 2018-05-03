#!/bin/sh

# ---------------
# Initial parameters
# ---------------

QEDIR="mpirun -np 4 /qe-path/bin"
OUT_DIR="./"

# ---------------
# Run
# ---------------

echo "SCF"
$QEDIR/pw.x < MgO.scf.in > $OUT_DIR/MgO.scf.out 2> /dev/null

echo "BAND"
$QEDIR/pw.x < MgO.band.in > $OUT_DIR/MgO.band.out 2> /dev/null

echo "BANDS"
$QEDIR/bands.x < MgO.bands.in > $OUT_DIR/MgO.bands.out 2> /dev/null

echo "NSCF"
$QEDIR/pw.x < MgO.nscf.in > $OUT_DIR/MgO.nscf.out 2> /dev/null

echo "DOS"
$QEDIR/dos.x < MgO.dos.in > $OUT_DIR/MgO.dos.out 2> /dev/null




