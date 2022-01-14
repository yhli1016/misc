#! /bin/bash

# Last modifiled on 2016-07-11
#
# This program extracts the stress at each ionic step during vc-relax calculation.
# The stress information is useful when checking the convergence and Pulay strees.
#
# Usage: stress.sh vc-relax.out

outfile=$1

line_number_list=$(cat -n $outfile | grep 'total   stress' | awk '{print $1}')
for line_number in $line_number_list; do
	awk '(NR>='$line_number'+1)&&(NR<='$line_number'+3) \
	    {printf "%10.2f%10.2f%10.2f\n", $4, $5, $6}' $outfile
	echo
done
