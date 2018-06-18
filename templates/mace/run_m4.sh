#! /bin/bash

topdir=$(pwd)

cd $topdir/01-scf
m4 scf.m4 -I../share | awk 'NF>0' > scf.in

for pref in 02-wfn 03-wfnq 04-wfn_path 05-wfn_fi; do
	cd $topdir/$pref
	m4 bands.m4 -I../share | awk 'NF>0' > bands.in
	m4 p2b.m4 -I../share | awk 'NF>0' > p2b.in
done
