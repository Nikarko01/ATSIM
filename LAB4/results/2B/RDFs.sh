#!/bin/bash -f

# Parameters
SCRIPTPATH=../../scripts/rdf.py
DATAPATH=../2A

# Start loop
for json in $DATAPATH/*.json
do
jout=$(basename $json) # exctracting file name
jout=${jout%.json} # remove extension
echo "Computing Radial Distribution Function (RDF) for "$json"..."
python $SCRIPTPATH $json -o $jout.dat --init-time 2.0 --n-bins 100 --stepsize-t 10 > slurm.out &
echo $jout".dat written"
done

echo "RDF finished!"
echo "Log output to slurm.out"
