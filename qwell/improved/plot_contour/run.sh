#! /bin/bash

date > eps1.dat
mpirun -np 4 python gen_eps1.py >> eps1.dat
wait

date > eps2.dat
mpirun -np 4 python gen_eps2.py >> eps2.dat
wait
