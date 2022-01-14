#! /bin/bash

# This program extracts the maximum for given band from eqp1.dat produced by sigma.$flavor.x
#
# This program is useful when calling plotband.sh.
#
# Usage: ef_bgw.sh eqp1.dat bndindex

bndfile=$1
bndindex=$2

ef_mf=$(awk '($1==1)&&($2=='$bndindex') {print $3}' $bndfile | sort -g | tail -1)
ef_gw=$(awk '($1==1)&&($2=='$bndindex') {print $4}' $bndfile | sort -g | tail -1)

printf "%16.9f%16.9f\n" $ef_mf $ef_gw
