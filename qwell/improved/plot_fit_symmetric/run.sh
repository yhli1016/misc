#! /bin/bash

date >> theo_co.dat
mpirun -np 4 python gen_theo_co.py >> theo_co.dat
wait
date >> theo_co.dat

date >> theo_fi.dat
mpirun -np 4 python gen_theo_fi.py >> theo_fi.dat
wait
date >> theo_fi.dat
