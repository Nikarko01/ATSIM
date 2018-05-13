#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 24
#SBATCH --cpus-per-task 1
#SBATCH --time=00:30:00
#SBATCH --mem=14000
# SBATCH --partition=mse-468
#SBATCH --account=mse-468
# SBATCH --reservation=mse-468

#These are the libraries necessary for our code to run

##### QE-Path ######
PW_LAUNCH="/home/nbarchi/Codes/espresso-5.1-hubbard-master_28.02.2018/bin/pw.x"

PSEUDODIR="../PP"
OUTODIR="../tmp"
PSEUDO='C.revpbe-rrkjus.adc.UPF'
# Input data
#
# List of values for volumes (The cell contains 2 atoms!)
LISTVOL="293.3 292.6 292.0 291.4 290.7 290.1 289.4 288.8 288.1 287.5" # TO COMPLETE

# List of values for the c/a ratio 
LISTC_O_A="3.094 3.088 3.081 3.074 3.067 3.060 3.054 3.047 3.040 3.033" # TO COMPLETE

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

echo "$A $CoA $vol $EN" >> $PREFIX"_E_vs_V_acoa.dat"

echo "Energy $EN, Vol =$vol"
done
done
