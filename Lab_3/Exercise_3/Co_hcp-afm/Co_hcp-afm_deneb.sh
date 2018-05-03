#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 8
#SBATCH --cpus-per-task 1
#SBATCH --time=00:30:00
#SBATCH --mem=14000
#SBATCH --account=mse-468
#SBATCH --reservation=mse-468-04-10

#These
module purge
module load intel/17.0.2
module load intel-mpi/2017.2.174
module load intel-mkl/2017.2.174
module load fftw/3.3.6-pl2
module load espresso/6.1.0-mpi

PW_LAUNCH="srun pw.x"

# Input data
#
# List of values for the lattice parameter
LISTA=" 4.65 " # TO COMPLETE
# List of values for the c/a ratio 
LISTC=" 1.62 " # TO COMPLETE

# Prefix name
PREFIX="Co_hcp-afm"

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
    outdir = '../temp/'
    pseudo_dir = '../PP/'
 /        
 &system    
    ibrav=  4
    celldm(1) = $a
    celldm(3) = $c_over_a
    nat=  2
    ntyp= 2
    ecutwfc = 55
    ecutrho = 440
    occupations = 'smearing'
    degauss = 0.01
    smearing = 'm-v'
    nspin = 2
    starting_magnetization(1) =  0.7
    starting_magnetization(2) =  -0.7
 /
 &electrons
    mixing_beta = 0.7 
    conv_thr =  1.0d-8
 /
 ATOMIC_SPECIES
  CoU  58.933194  Co.pbe-n-rrkjus_psl.1.0.0.UPF
  CoD  58.933194  Co.pbe-n-rrkjus_psl.1.0.0.UPF
 ATOMIC_POSITIONS crystal
  CoU  0.333333 0.666667 0.25
  CoD  0.666667 0.333333 0.75
 K_POINTS automatic
    12 12 6 0 0 0 
EOF

$PW_LAUNCH < $INFILE > $OUTFILE

EN=` grep ! $OUTFILE | egrep -o "([+-])?[0-9]+(\.[0-9]+)?" `

echo "$c_over_a $a $EN" >> ${PREFIX}_results.dat

echo "c/a=${c_over_a}, a=${a} Bohr, Energy=${EN} Ry"

done
done
