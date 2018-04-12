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


LISTA="7.97"             	# List of values of lattice parameter to try
LISTECUT="xxx"              # List of plane-wave cutoffs to try
LISTK="1 2 3 4 5 6 7 8"  	# List of number of k-points per dimension to try.


# Files of interest:
TMP_DIR="./tmp"         	# where temporary data will be stored.
PSEUDO_DIR="../pseudo"  	# where pseudopotentials are stored.
OUT_DIR="./results"     	# where input and output will be
                        	# created once the script runs.


# check whether ECHO has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi


# Security checks:
if [ ! -d $TMP_DIR ]; then
   mkdir $TMP_DIR
fi
if [ ! -d $OUT_DIR ]; then
   mkdir $OUT_DIR
fi


# Output header
$ECHO "Energy\tForce\tPressure\tEcutwf\tKmesh" >> data

# Start loops on plane-wave cutoffs, k-point grids, and lattice constants:
for ecut in $LISTECUT 
do 
for k in $LISTK 
do 
for a in $LISTA 
do

#### Choose on how many processors to run depending on the system size

PW_LAUNCH="srun pw.x"

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
         Mg 0.05 0.05 0.05
         O  0.50 0.50 0.50
      K_POINTS {automatic}
         $k $k $k  0 0 0
EOF

# Run PWscf to create new output file:
$ECHO " running the scf calculation for..." MgO.scf.a=$a.ecut=$ecut.k=$k.in
$PW_LAUNCH < $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.in > $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.out

# Extract data
E=$(grep ! $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.out | awk '{print $5}')
F=$(grep "Total force =" $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.out | awk '{print $4}')
P=$(grep "P=" $OUT_DIR/MgO.scf.a=$a.ecut=$ecut.k=$k.out | awk '{print $6}')
$ECHO "$E\t$F\t$P\t$ecut\t$k" >> data

# Finish loops on plane-wave cutoffs, k-point grids, and lattice constants:
done
done
done

# Clean up:
rm -rf $TMP_DIR/*