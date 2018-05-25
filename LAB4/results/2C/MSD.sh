#!/bin/bash -f

# Parameters
SCRIPTPATH=../../../scripts/msd.py
DATAPATH=../../2A

# Start loop
for json in $DATAPATH/*.json
do
jout=$(basename $json) # exctracting file name
jout=${jout%.json} # remove extension
echo "Computing Radial Distribution Function (RDF) for "$json"..."
python $SCRIPTPATH $json -o $jout.dat --stepsize-t 6 > slurm.out &
echo $jout".dat written"
done

echo "RDF finished!"
echo "Log output to slurm.out"
