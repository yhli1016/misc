#! /bin/bash

date >> fmin.out
mpirun -np 4 python fmin.py >> fmin.out
wait
date >> fmin.out

date >> leastsq.out
mpirun -np 4 python leastsq.py >> leastsq.out
wait
date >> leastsq.out
