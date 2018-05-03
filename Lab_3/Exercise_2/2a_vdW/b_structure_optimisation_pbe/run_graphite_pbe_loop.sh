#!/bin/sh

# ---------------
# Initial parameters
# ---------------

PWROOT="mpirun -np 1 YOUR_PATH_TO_pw.x"
PSEUDODIR="../PP"
OUTODIR="../tmp"

# ---------------
# Run parameters
# ---------------

THISNAME='graphite_eos_pbe'
PSEUDO='C.pbe-rrkjus.UPF'

ALIST=""

CLIST=""

for A in $ALIST
do
for C in $CLIST
do

CoA=`echo "scale=8;$C/$A" | bc`
VAL1=`echo "scale=8;(1/4)*$CoA" | bc`
VAL2=`echo "scale=8;(3/4)*$CoA" | bc`

echo "Running $THISNAME.a.$A.c.$C: a = $A, c=$C..."

cat <<EOF > $THISNAME.a.$A.c.$C.scf.in 
&control
  calculation  = 'scf'
  restart_mode = 'from_scratch'
  prefix       = '$THISNAME'
  pseudo_dir   = '$PSEUDODIR'
  outdir       = '$OUTODIR'
  wf_collect   = .true.
 /
&system
  nat=4, ntyp=1,
  ibrav=4, celldm(1)=$A, celldm(3)=$CoA
  ecutwfc = 40.0
  ecutrho = 320.0
/
&electrons
  diagonalization = 'cg'
  conv_thr = 1.0e-8
/
ATOMIC_SPECIES
 C 12.0107 $PSEUDO

ATOMIC_POSITIONS (crystal)
 C  0.33333333333333  0.66666666666666  0.25000000000000
 C  0.66666666666666  0.33333333333333  0.25000000000000
 C  0.33333333333333  0.66666666666666  0.75000000000000
 C  0.00000000000000  0.00000000000000  0.75000000000000

K_POINTS {automatic}
4 4 2 1 1 1
EOF

$PWROOT/bin/pw.x < $THISNAME.a.$A.c.$C.scf.in > $THISNAME.a.$A.c.$C.scf.out 2> /dev/null

EN=`cat $THISNAME.a.$A.c.$C.scf.out | grep -e ! | egrep -o "([+-])?[0-9]+(\.[0-9]+)?"`

echo "$A $C $EN"   >> $THISNAME"_ac.dat"
echo "$A $CoA $EN" >> $THISNAME"_acoa.dat"

echo "Ok - Energy $EN"
done
done

