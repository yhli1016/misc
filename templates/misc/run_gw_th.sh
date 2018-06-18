#! /bin/bash

# mandatory for BGW
source /opt/intel/composer_xe_2013.0.079/bin/iccvars.sh intel64
source /opt/intel/composer_xe_2013.0.079/bin/ifortvars.sh intel64

ncore=28
queue=TH_HPC1
MPIRUN="yhrun -n $ncore -p $queue"
topdir=$(pwd)
