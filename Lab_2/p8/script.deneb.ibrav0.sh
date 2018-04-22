#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 8
#SBATCH --time=00:30:00

#These
module purge
module load intel/17.0.2
module load intel-mpi/2017.2.174
module load intel-mkl/2017.2.174
module load fftw/3.3.6-pl2
module load espresso/6.1.0-mpi


LISTX="0.0 0.000714 0.001429 0.002143 0.002857 0.003571 0.004286 0.005 0.005714 0.006429 0.007143 0.007857 0.008571 0.009286 0.01" # List of values of lattice parameter to try
LISTECUT="70"          # List of plane-wave cutoffs to try
LISTK="4"               # List of number of k-points per dimension to try.


# Files of interest:
TMP_DIR="/scratch/nsbarchi/tmp"         # where temporary data will be stored.
PSEUDO_DIR="../pseudo"  # where pseudopotentials are stored.
OUT_DIR="./results.8B.ibrav0"     # where input and output will be
                        # created once the script runs.

#------------------------------------------------------------------------------
# Checks
#------------------------------------------------------------------------------

# check whether ECHO has the -e option --SECURITY_CHECK_1
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

# Security checks: create a new $TMP_DIR --SECURITY_CHECK_2
count=0
testdir=$TMP_DIR
while [ -d "$testdir" ]
do
   $ECHO "$testdir already exists"
   let "count++"
   testdir=$TMP_DIR
   testdir+=$($ECHO "_"$count)
done
$ECHO "Found available name: $testdir"
TMP_DIR=$testdir

# Make dir $OUT_DIR --SECURITY_CHECK_3
if [ ! -d $OUT_DIR ]; then
   mkdir $OUT_DIR
fi


# Output header
$ECHO "Energy\t Ecutwf \t Deform" >> $OUT_DIR/data

# Start loops on plane-wave cutoffs, k-point grids, and lattice constants:
for ecut in $LISTECUT 
do 
for k in $LISTK 
do 
for x in $LISTX 
do

#### Choose on how many processors to run depending on the system size

PW_LAUNCH="srun pw.x"
alat=8.04475
a=$(echo "1+$x"            | bc -l)
b=$(echo "1-$x"            | bc -l)
c=$(echo "1+$x^2/(1-$x^2)" | bc -l)

# Create new input file:
cat > $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.in << EOF
      &CONTROL
         calculation = 'relax'
         restart_mode = 'from_scratch'
         prefix = 'diamond.$a.$ecut.$k'
         tstress = .true.
         tprnfor = .true.
         pseudo_dir = '$PSEUDO_DIR'
         outdir='$TMP_DIR'
         etot_conv_thr = 1.0D-5 ! 1 convergence criteria  
         forc_conv_thr = 1.0D-4 ! 2nd conv. criteria
      /

      &SYSTEM
         ibrav = 0
         celldm(1) = $alat



         nat = 8
         ntyp = 2
         ecutwfc = $ecut
      /

      &ELECTRONS
         diagonalization = 'david'
         mixing_mode = 'plain'
         mixing_beta = 0.7
         conv_thr = 1.0d-8
      /
      
      &IONS
         ion_dynamics = 'bfgs'
      /

      &CELL
         cell_dynamics = 'bfgs'
         press = 0.0d0
         press_conv_thr = 0.01D0
      /

      ATOMIC_SPECIES
         Mg  24.305   Mg.pbe.UPF
         O   15.9994  O.pbe.UPF
      CELL_PARAMETERS (alat= 8.04475)
             $a        0.000000000   0.000000000
         0.000000000       $b        0.000000000
         0.000000000   0.000000000       $c     
      ATOMIC_POSITIONS {alat} 
         Mg 0.00 0.00 0.00
         Mg 0.00 0.50 0.50
         Mg 0.50 0.00 0.50
         Mg 0.50 0.50 0.00
         O  0.50 0.00 0.00
         O  0.00 0.50 0.00
         O  0.00 0.00 0.50
         O  0.50 0.50 0.50
      K_POINTS {automatic}
         $k $k $k  0 0 0
EOF

# Run PWscf to create new output file:
$ECHO " running the scf calculation for..." MgO.scf.a=$a.ecut=$ecut.k=$k.in
$PW_LAUNCH < $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.in > $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.out

# Extract data
E=$(grep "Final energy" $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.out | awk '{print $4}')
$ECHO "$E\t$ecut\t$x" >> $OUT_DIR/data

# Finish loops on plane-wave cutoffs, k-point grids, and lattice constants:
done
done
done

# Clean up:
rm -rf $TMP_DIR/*
