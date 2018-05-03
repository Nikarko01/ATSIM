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

# ---------------
# Initial parameters
# ---------------

SRUN="srun "
OUT_DIR="./"

# ---------------
# Run
# ---------------

echo "SCF"
$SRUN pw.x < MgO.scf.in > $OUT_DIR/MgO.scf.out

echo "BAND"
$SRUN pw.x < MgO.band.in > $OUT_DIR/MgO.band.out 

echo "BANDS"
$SRUN bands.x < MgO.bands.in > $OUT_DIR/MgO.bands.out 

echo "NSCF"
$SRUN pw.x < MgO.nscf.in > $OUT_DIR/MgO.nscf.out 

echo "DOS"
$SRUN dos.x < MgO.dos.in > $OUT_DIR/MgO.dos.out 




