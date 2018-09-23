#! /bin/bash

date >> fmin.out
mpirun -np 5 python fmin.py >> fmin.out
wait
date >> fmin.out

date >> leastsq.out
mpirun -np 5 python leastsq.py >> leastsq.out
wait
date >> leastsq.out
