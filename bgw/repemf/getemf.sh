#! /bin/bash

# Sometimes GW and BSE calculations are performed on different hosts. In this
# case, two WFN_co files will be created. However, small numerical errors are
# unavoidable due to the difference in hardware and software configurations of
# different hosts, which will cause inteqp and absorption codes refuse to work.
#
# This program solves the problem by replacing the meanfield energies in eqp1.dat
# obtained on another host with the accurate energies in 02-wfn/xxx.save/K0000X/
# eigenval.xml obtained on the current host.
#
# This script is an auxiliary program that extracts meanfield energies. Invoke it
# before calling repemf.x.
#
# Usage: getemf.sh ikmax, where ikmax is the maximum number of kpoints.

ikmax=$1

if [ -f emf.dat ]; then
	rm emf.dat
fi
touch emf.dat

for ik in $(seq $ikmax); do
	pref=$(printf "K%05d" $ik)
	nl0=$(awk '/EIGENVAL/ {print NR}' $pref/eigenval.xml | head -1)
	nl1=$(awk '/EIGENVAL/ {print NR}' $pref/eigenval.xml | tail -1)
	awk '(NR>'$nl0') && (NR<'$nl1') {print}' $pref/eigenval.xml >> emf.dat
done
