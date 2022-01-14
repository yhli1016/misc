#! /bin/bash

# NAMD assumes that directories under 'rundir' should start from 0001,
# but we may start from the 1001st frame of the ordinary MD trajectory.
# So we need to rename the directories before running NAMD.
#
# This program should be run in the directory defined by 'rundir'.

nmin=1001
nmax=2000
dn=1000

for i in $(seq $nmin $nmax); do
        dir1=$(echo $i | awk '{printf "%04d", $1}')
        dir2=$(echo $i $dn | awk '{printf "%04d", $1 - $2}')
        echo "Moving" $dir1 "to" $dir2
        mv $dir1 $dir2
done
