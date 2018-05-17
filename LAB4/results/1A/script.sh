#!/bin/bash -f

# Script written to run GULP multiple times on one node.
# The script can loop over 4 different input parameters:
### the time step, specified in list_time_steps
### The temperature, specified in list_temperatures
### The supercell size, specified in list_supercell_size
### The equilibration time, specified in list_equilibration_time

# There are several ways to set your list in BASH.
### For explicit definition of e.g. your timesteps, you can do
##### list_time_steps="0.001 0.003 0.01" #######

### You can also you the seq command to create a sequence as in:
##### list_time_steps=`seq 0.001 0.001 0.01`
### This creates a sequence of values between 0.001 and 0.01, every 0.001

list_time_steps='0.05 0.1 0.5 0.7 1 2 5 10'
list_temperatures='3000'
list_supercell_size='2'

# This is the executable path for gulp
exec="/home/nbarchi/Codes/gulp-5.0/Src/gulp"

# This is the executable path for the parser
parser="python /home/nbarchi/ATSIM/LAB4/scripts/parser.py"

# These are constants that the script does not loop over, for production time
# and sampling stepsize
prod="10"
sample="0.02"
equilibration="5.0"
# Start loops
for cellsize in $list_supercell_size; do
    for temperature in $list_temperatures; do
            for time_step in $list_time_steps; do
                echo "Cell size ${cellsize} Temperature ${temperature}, Equilibration time ${equilibration}, Timestep $time_step ps, Production $prod ps, Sampling $sample"
                base_name="${cellsize}_${temperature}_${time_step}"
                cat > md_test_${base_name}.gin << EOF
conv md
title
molecular dynamics of silver
end
cell
    4.085    4.085     4.085    90 90 90
fractional    4
Ag     0.000000000     0.000000000     0.000000000
Ag     0.000000000     0.500000000     0.500000000
Ag     0.500000000     0.000000000     0.500000000
Ag     0.500000000     0.500000000     0.000000000
#
#  modified version of the Cleri-Rosato potential for GULP (Ag)
#  from F. Cleri and V. Rosato
#  Phys. Rev. B, 48, 22 (1993)
#
#  This potential uses the Embedded Atom Model
#
species
Ag core  0.000
manybody
Ag core Ag core  0.0 8.0
#
#  Functional
#
eam_functional square_root
Ag core  1.000
#
#  Density term
#
eam_density exponential 0
Ag core  1.387684 2.173423 2.888531
#
#  Repulsive two-body components -
#  second power is a dummy argument
#
buckingham
Ag core Ag core    11454.950282 0.264324 0.0 0.0 12.0
#
# MD specific parameters
#
supercell ${cellsize} ${cellsize} ${cellsize}
ensemble nve
temperature ${temperature}
integrator leapfrog verlet
equilbration ${equilibration} ps
production   ${prod} ps
timestep     ${time_step} ps
sample       ${sample} ps
write        ${sample} ps
output trajectory ascii Trajectories_${base_name}.trg
dump every 100 Recovery_${base_name}.res
EOF

            $exec < md_test_${base_name}.gin > md_test_${base_name}.gout
            $parser Trajectories_${base_name}.trg md_test_${base_name}.gout output_${base_name} -e E_of_t_${base_name}.dat

        done
    done
done
