#! /bin/bash

topdir=$(pwd)
launcher="mpirun -np 20"

cd $topdir/01-scf
$launcher pw.x -in scf.in &> scf.out
wait

for pref in 02-wfn 03-wfnq 04-wfn_path; do
    cd $topdir/$pref
    cp -r ../01-scf/*.save .
    wait
    $launcher pw.x -in bands.in &> bands.out
    wait
    $launcher pw2bgw.x -in p2b.in &> p2b.out
    wait
    rm -rf *.save *wfc* *igk*
    wait
done
