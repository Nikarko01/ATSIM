#!/bin/bash

# Input data:
LISTA="7.969"            # List of values of lattice parameter to try
LISTECUT="20 25"        # List of plane-wave cutoffs to try
LISTK="2 4"             # List of number of k-points per dimension to try

# Files of interest:
TMP_DIR="./tmp"        # where temporary data will be stored
PSEUDO_DIR="./pseudo"  # where pseudopotentials are stored
OUT_DIR="./Test_script"       # where input and output will be
                            # created once the script runs.

PW_LAUNCH='pw.x'            # This is QE executable file.

# Security checks:
if [ ! -d $TMP_DIR ]; then
   mkdir $TMP_DIR
fi

if [ ! -d $OUT_DIR ]; then
   mkdir $OUT_DIR
fi

# Start loops on plane-wave cutoffs, k-point grids, and lattice constants:
for ecut in $LISTECUT
do
for k in $LISTK
do
for a in $LISTA
do

echo 'Doing ecut, k, a: ' $ecut $k $a
# Create new input file:
cat > $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.in << EOF
      &control
         calculation = 'scf'
         restart_mode = 'from_scratch'
         prefix = 'diamond.$a.$ecut.$k'
         tstress = .true.
         tprnfor = .true.
         pseudo_dir = '$PSEUDO_DIR'
         outdir='$TMP_DIR'
      /
      &system
         ibrav = 2
         celldm(1) = $a
         nat = 2
         ntyp = 2
         ecutwfc = $ecut
      /
      &electrons
         diagonalization = 'david'
         mixing_mode = 'plain'
         mixing_beta = 0.7
         conv_thr = 1.0d-8
      /
      ATOMIC_SPECIES
         Mg  24.305   Mg.pbe.UPF
         O   15.9994  O.pbe.UPF
      ATOMIC_POSITIONS {alat} 
         Mg 0.00 0.00 0.00
         O  0.50 0.50 0.50
      K_POINTS {automatic}
         $k $k $k  0 0 0
EOF

# Run PWscf to create new output file:
$PW_LAUNCH < $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.in > $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.out

# Finish loops on plane-wave cutoffs, k-point grids, and lattice constants:
done
done
done

rm -r $TMP_DIR/*
