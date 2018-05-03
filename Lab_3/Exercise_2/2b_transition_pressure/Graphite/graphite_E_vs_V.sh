#!/bin/bash -f


##### QE-Path ######
PW_LAUNCH="mpirun -np 16 YOUR_PATH_TO_pw.x "

PSEUDODIR="../PP"
OUTODIR="../tmp"
PSEUDO='C.revpbe-rrkjus.adc.UPF'
# Input data
#
# List of values for volumes (The cell contains 2 atoms!)
LISTVOL="" # TO COMPLETE
# List of values for the c/a ratio 
LISTC_O_A="" # TO COMPLETE

# Prefix name
PREFIX="graphite"

#Loop over volumes
for vol in $LISTVOL
do
# Loop over c/a ratios 
for CoA in $LISTC_O_A
do

# Recover the lattice parameter from the volume and the c/a ratio
A=`echo "$vol $CoA" | awk '{printf "%.7f", (2*$1/sqrt(3)/$2)^(1.0/3)}'`

INFILE="$PREFIX-scf-$vol-$CoA-$A.in"
OUTFILE="$PREFIX-scf-$vol-$CoA-$A.out"
rm -f $INFILE
rm -f $OUTFILE

cat > $INFILE << EOF
&control
  calculation  = 'scf'
  restart_mode = 'from_scratch'
  prefix       = '$PREFIX'
  pseudo_dir   = '$PSEUDODIR'
  outdir       = '$OUTODIR'
  wf_collect   = .true.
 /
&system
  nat=4, ntyp=1,
  ibrav=4, celldm(1)=$A, celldm(3)=$CoA
  ecutwfc = 40.0
  ecutrho = 320.0
  input_dft = 'VDW-DF'
/
&electrons
  diagonalization = 'cg'
  conv_thr = 1.0e-8
/
ATOMIC_SPECIES
 C 12.0107 $PSEUDO

ATOMIC_POSITIONS (crystal)
 C  0.33333333333333  0.66666666666666  0.25000000000000
 C  0.66666666666666  0.33333333333333  0.25000000000000
 C  0.33333333333333  0.66666666666666  0.75000000000000
 C  0.00000000000000  0.00000000000000  0.75000000000000

K_POINTS {automatic}
4 4 2 1 1 1
EOF

$PW_LAUNCH < $INFILE > $OUTFILE

EN=`cat $OUTFILE | grep -e ! | egrep -o "([+-])?[0-9]+(\.[0-9]+)?"`

echo "$A $CoA $vol $EN" >> $PREFIX"_acoa.dat"

echo "Energy $EN, Vol =$vol"
done
done
