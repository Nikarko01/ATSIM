#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --cpus-per-task 1
#SBATCH --time=00:30:00
# SBATCH --mem=14000
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

# ---------------
# Initial parameters
# ---------------

PWROOT="srun -n 24 pw.x"
PSEUDODIR="../PP"
OUTODIR="../tmp"

# ---------------
# Run parameters
# ---------------

THISNAME='graphite_eos_pbe'
PSEUDO='C.pbe-rrkjus.UPF'

ALIST="2.088 2.169 2.251 2.333 2.415 2.497 2.579 2.661 2.743 2.824"

CLIST="5.692 5.835 5.979 6.122 6.266 6.409 6.553 6.696 6.839 6.983 7.126 7.270 7.413 7.557 7.700"

for A in $ALIST
do
for C in $CLIST
do
rm -rf $OUTODIR/*
CoA=`echo "scale=8;$C/$A" | bc`
VAL1=`echo "scale=8;(1/4)*$CoA" | bc`
VAL2=`echo "scale=8;(3/4)*$CoA" | bc`

echo "Running $THISNAME.a.$A.c.$C: a = $A, c=$C..."

cat <<EOF > $THISNAME.a.$A.c.$C.scf.in 
&control
  calculation  = 'scf'
  restart_mode = 'from_scratch'
  prefix       = '$THISNAME'
  pseudo_dir   = '$PSEUDODIR'
  outdir       = '$OUTODIR'
  wf_collect   = .true.
 /
&system
  nat=4, ntyp=1,
  ibrav=4, celldm(1)=$A, celldm(3)=$CoA
  ecutwfc = 40.0
  ecutrho = 320.0
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

$PWROOT < $THISNAME.a.$A.c.$C.scf.in > $THISNAME.a.$A.c.$C.scf.out 2> /dev/null

EN=`cat $THISNAME.a.$A.c.$C.scf.out | grep -e ! | egrep -o "([+-])?[0-9]+(\.[0-9]+)?"`

echo "$A $C $EN"   >> $THISNAME"_ac.dat"
echo "$A $CoA $EN" >> $THISNAME"_acoa.dat"

echo "Ok - Energy $EN"
done
done

