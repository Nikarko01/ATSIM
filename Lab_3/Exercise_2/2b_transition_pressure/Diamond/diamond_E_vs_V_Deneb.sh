#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 24
#SBATCH --cpus-per-task 1
#SBATCH --time=00:30:00
#SBATCH --mem=14000
#SBATCH --partition=mse-468
#SBATCH --account=mse-468
#SBATCH --reservation=mse-468

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
LISTVOL="" # TO COMPLETE

# Prefix name
PREFIX="diamond"

#Loop over volumes
for vol in $LISTVOL
do

# Recover the lattice parameter from the volume and the c/a ratio
A=`echo "$vol" | awk '{printf "%.7f", (4*$1)^(1.0/3)}'`

INFILE="$PREFIX-scf-$vol-$A.in"
OUTFILE="$PREFIX-scf-$vol-$A.out"
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
  nat=2, ntyp=1,
  ibrav=2, celldm(1)=$A
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

ATOMIC_POSITIONS (alat)
 C   0.0   0.0   0.0
 C   0.25  0.25  0.25
 
K_POINTS {automatic}
 8 8 8 1 1 1
EOF

$PW_LAUNCH < $INFILE > $OUTFILE

EN=`cat $OUTFILE | grep -e ! | egrep -o "([+-])?[0-9]+(\.[0-9]+)?"`

echo "$A $vol $EN" >> $PREFIX"_E_vs_V.dat"

echo "Energy $EN, Vol =$vol"
done