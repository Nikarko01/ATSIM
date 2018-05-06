#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=6:00:00

#source /ssoft/spack/bin/slmodules.sh -r stable -v
#module load intel/17.0.2
#module load fftw/3.3.6-pl2
#module load intel-mkl/2017.2.174
#module load intel-mpi/2017.2.174
#module load elpa/2016.05.004
#module load espresso/6.1.0-mpi
module load intel intel-mpi intel-mkl espresso


PW_LAUNCH="srun pw.x"

# Input data
#
# List of values for the lattice parameter
LISTA=" 6.035 6.04 6.045 6.05 6.045 6.04 6.035 6.035 " # TO COMPLETE
# List of values for the c/a ratio 
LISTC=" 1.59 1.60 1.61 1.615 1.62 1.625 1.63 1.635 1.64 " # TO COMPLETE

# Prefix name
PREFIX="Mg"

# Loop over a values
for a in $LISTA
do
# Loop over c/a ratios 
for c_over_a in $LISTC
do

INFILE="${PREFIX}_scf_${a}_${c_over_a}.in"
OUTFILE="${PREFIX}_scf_${a}_${c_over_a}.out"
rm -f $INFILE
rm -f $OUTFILE

cat > $INFILE << EOF
 &control
    calculation = 'scf'
    restart_mode='from_scratch'
    prefix='$PREFIX'
    tstress = .true.
    tprnfor = .true.
    outdir = './tmp'
    pseudo_dir = '/home/jaffrelo/PSEUDOS_GOOD/'
 /        
 &system    
    ibrav = 4
    celldm(1) = $a
    celldm(3) = $c_over_a 
    nat = 2
    ntyp = 1
    ecutwfc = 45.0
    ecutrho = 360.0
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.02
 /
 &electrons
    diagonalization = 'david'
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr = 1.0d-8
 /
 ATOMIC_SPECIES
  Mg  24.305   Mg.pbe.uspp.UPF
 ATOMIC_POSITIONS crystal
  Mg  0.333333 0.666667 0.25
  Mg  0.666667 0.333333 0.75
 K_POINTS automatic
    12 12 6 0 0 0 
EOF

$PW_LAUNCH < $INFILE > $OUTFILE

EN=` grep ! $OUTFILE | egrep -o "([+-])?[0-9]+(\.[0-9]+)?" `

echo "$c_over_a $a $EN" >> ${PREFIX}_results.dat

echo "c/a=${c_over_a}, a=${a} Bohr, Energy=${EN} Ry"

done
done

