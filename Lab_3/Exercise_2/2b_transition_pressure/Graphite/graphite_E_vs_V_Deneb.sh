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
module purge
module load intel/17.0.2
module load intel-mpi/2017.2.174
module load intel-mkl/2017.2.174
module load fftw/3.3.6-pl2
module load espresso/6.1.0-mpi


##### QE-Path ######
PW_LAUNCH="srun -n 24 pw.x"

PSEUDODIR="../PP"
OUTODIR="../tmp"
PSEUDO='C.revpbe-rrkjus.adc.UPF'
# Input data
#
# List of values for volumes (The cell contains 2 atoms!)
LISTVOL="12.82 12.79 12.76 12.73 12.71 12.68 12.65 12.62 12.59 12.56" # TO COMPLETE

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

echo "$A $CoA $vol $EN" >> $PREFIX"_acoa.dat"

echo "Energy $EN, Vol =$vol"
done
done
