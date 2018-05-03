#!/bin/bash -f

PW_LAUNCH="mpirun -np 1 /qe-path/pw.x"

# Input data
#
# List of values for the lattice parameter
LISTA=" 6.66 " # TO COMPLETE

# Prefix name
PREFIX="Co_fcc-fm"

# Loop over a values
for a in $LISTA
do

INFILE="$PREFIX.scf.$a.in"
OUTFILE="$PREFIX.scf.$a.out"
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
    ibrav= 2
    celldm(1) = $a
    nat=  1
    ntyp= 1
    ecutwfc = 55
    ecutrho = 440
    occupations = 'smearing'
    degauss = 0.01
    smearing = 'm-v'
    nspin = 2
    starting_magnetization(1) =  0.7
 /
 &electrons
    mixing_beta = 0.7 
    conv_thr =  1.0d-8
 /
 ATOMIC_SPECIES
  Co  58.933194  Co.pbe-n-rrkjus_psl.1.0.0.UPF 
 ATOMIC_POSITIONS crystal
  Co  0.0 0.0 0.0
 K_POINTS automatic
    12 12 12 0 0 0 
EOF

$PW_LAUNCH < $INFILE > $OUTFILE

EN=` grep ! $OUTFILE | egrep -o "([+-])?[0-9]+(\.[0-9]+)?" `

echo "$a $EN" >> ${PREFIX}_results.dat

echo "a=${a} Bohr, Energy=${EN} Ry"

done
