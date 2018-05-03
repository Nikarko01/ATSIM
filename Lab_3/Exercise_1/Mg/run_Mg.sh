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
$QEDIR/pw.x < Mg.scf.in > $OUT_DIR/Mg.scf.out 2> /dev/null

echo "BAND"
$QEDIR/pw.x < Mg.band.in > $OUT_DIR/Mg.band.out 2> /dev/null

echo "BANDS"
$QEDIR/bands.x < Mg.bands.in > $OUT_DIR/Mg.bands.out 2> /dev/null

echo "NSCF"
$QEDIR/pw.x < Mg.nscf.in > $OUT_DIR/Mg.nscf.out 2> /dev/null

echo "DOS"
$QEDIR/dos.x < Mg.dos.in > $OUT_DIR/Mg.dos.out 2> /dev/null




