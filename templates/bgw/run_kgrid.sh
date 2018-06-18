#! /bin/bash

# 02-wfn
nk1=
nk2=
nk3=

dk1=0
dk2=0
dk3=0

dq1=0
dq2=0
dq3=0

# 03-wfn
dq1q="0.001"
dq2q=0
dq3q=0

# 05-wfn_fi
nk1_fi=
nk2_fi=
nk3_fi=

# fft grid
ng1=
ng2=
ng3=

# pd file
pdfile=

#-----------------------------------------------------------------------------------------------------
pw2kgrid.py $pdfile wfn.inp    $nk1    $nk2    $nk3    $dk1 $dk2 $dk3 $dq1  $dq2  $dq3  $ng1 $ng2 $ng3
pw2kgrid.py $pdfile wfnq.inp   $nk1    $nk2    $nk3    $dk1 $dk2 $dk3 $dq1q $dq2q $dq3q $ng1 $ng2 $ng3
pw2kgrid.py $pdfile wfn_fi.inp $nk1_fi $nk2_fi $nk3_fi $dk1 $dk2 $dk3 $dq1  $dq2  $dq3  $ng1 $ng2 $ng3

kgrid.x wfn.inp    wfn.out    wfn.log
kgrid.x wfnq.inp   wfnq.out   wfnq.log
kgrid.x wfn_fi.inp wfn_fi.out wfn_fi.log
