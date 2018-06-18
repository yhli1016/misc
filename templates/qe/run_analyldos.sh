#! /bin/bash

function analyldos()
{
	if [ ! -e trim.dat ]; then
		nl0=$(awk '(NF==2)&&($1=="F")&&($2=="F") {print NR -1}' proj.dat)
		awk 'NR>='$nl0' {print}' proj.dat > trim.dat
	fi

	cat > o.inp << EOF
trim.dat
1 161
1 60
$1
$(seq 5 $2)
EOF
	analyldos.x o.inp o.out
	wait
	analyldos_pp.sh o.out bndstr.pd
	wait
}

topdir=$(pwd)

for i in $(seq 1 10); do

	cd $topdir/$i

	idomax=$(echo $i | awk '{print $1 + 4}')

	analyldos $i $idomax

done
