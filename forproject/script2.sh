#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=12:00:00

source /ssoft/spack/bin/slmodules.sh -r stable -v
module load intel/17.0.2
module load fftw/3.3.6-pl2
module load intel-mkl/2017.2.174
module load intel-mpi/2017.2.174
module load elpa/2016.05.004
module load espresso/6.1.0-mpi

cat > scfcalc.in << EOF

 &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='BiFeO3_smearing'
    pseudo_dir = '/home/jaffrelo/PSEUDOS_GOOD/'
    outdir='./tmp'
    verbosity='high'
    tprnfor = .true.
    tstress = .true.
 /
 &system
    ibrav = 5
    A = 5.55791
    B = 5.55791
    C = 5.55791
    cosAB = 0.5026
    cosAC = 0.5026
    cosBC = 0.5026
    nspin = 2 
    nat = 10 
    ntyp = 4
    ecutwfc = 80
    ecutrho = 640
    occupations = 'smearing'
    tot_magnetization = 0
    starting_magnetization(1) = 0.5
    starting_magnetization(2) = -0.5
    smearing = 'mv'
    degauss = 0.01
 /
 &electrons
    conv_thr =  1.d-8
    mixing_beta = 0.7
 /
ATOMIC_SPECIES
Fe1  56.0   Fe.pbesol-spn-kjpaw_psl.0.2.1.UPF
Fe2  56.0   Fe.pbesol-spn-kjpaw_psl.0.2.1.UPF
Bi   83.0   bi_pbesol_v1.uspp.F.UPF
O    16.0   O.pbesol-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS {crystal}
 Fe1  0.00283   0.00283   0.00283
 Fe2  0.50284   0.50284   0.50283
 Bi   0.27560   0.27560   0.27560
 Bi   0.77561   0.77561   0.77561
 O    0.66880   0.21582   0.81432
 O    0.21582   0.81431   0.66880
 O    0.81432   0.66880   0.21582
 O    0.31434   0.71580   0.16880
 O    0.16880   0.31434   0.71580
 O    0.71581   0.16880   0.31434

K_POINTS {automatic}
 6 6 6  0 0 0

EOF

srun pw.x < scfcalc.in > scfcalc.out
