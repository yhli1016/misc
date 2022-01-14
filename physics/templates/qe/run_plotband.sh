#! /bin/bash

function plotband()
{
	cat > plotband.in << EOF
bndstr.dat
bndstr.pd
1
$1
1.889726125
EOF
	plotband.x < plotband.in > plotband.out
	wait
}

topdir=$(pwd)

for i in $(seq 1 10); do

	cd $topdir/$i

	vbm=$(grep occupied scf.out | awk '{print $5}')

	plotband $vbm
done
