#!/bin/bash


# Parameter
dir=./p5
# check whether ECHO has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

# Set function
min_number() {
	printf "%s\n" "$@" | sort -g | head -n1
}

$ECHO "Energy\t\tEcutwf\tAlat" >> $dir/data
$ECHO "#" $(date) >> $dir/data
for filename in $dir/results/*.in
do
	# ----------------------------------------------------------------------------
	# –––––––––––––––––– Extrapolate name for *in *out ---------------------------
	name=${filename%.*}
	$ECHO $filename
    # ----------------------------------------------------------------------------
	# –––––––––––––––––– Extrapolate variables -----------------------------------
	ec=$(grep ecutwfc $filename            | awk '{print $3}')
	# ep=$(grep ecutrho $filename            | awk '{print $3}')
	# k=$(grep -A 1 K_POINTS $filename | tail -n 1 | awk '{print $1" "$2" "$3}')
	# k=$(min_number $k)
	# sigma=$(grep degauss $filename            | awk '{print $3}')
	E=$(grep ! $name.out | awk '{print $5}')
	# F=$(grep "Total force =" $name.out | awk '{print $4}')
	# P=$(grep "P=" $name.out | awk '{print $6 $7}' | cut -d '=' -f2)
	a=$(grep 'celldm(1)' $filename            | awk '{print $3}')
	$ECHO "$E\t$ec\t$a" >> $dir/data
done