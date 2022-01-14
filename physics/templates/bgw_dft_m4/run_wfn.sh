#! /bin/bash

topdir=$(pwd)
RUN="mpirun -np $NP"

cd $topdir/01-scf
$RUN pw.x -nk 2 -in scf.in &> scf.out
wait

for pref in 02-wfn 03-wfnq 05-wfn_fi; do
    cd $topdir/$pref
    cp -r ../01-scf/*.save .
    wait
    $RUN pw.x -nk 2 -in bands.in &> bands.out
    wait
    $RUN pw2bgw.x -in p2b.in &> p2b.out
    wait
    rm -rf *.save *wfc* *igk*
    wait
done
