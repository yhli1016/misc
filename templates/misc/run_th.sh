#! /bin/bash

ncore=28
queue=TH_HPC1
MPIRUN="yhrun -n $ncore -p $queue"
topdir=$(pwd)
