#! /bin/bash

# Procedure of generating wave functions for NAMD calculation
#
# 1. Create the top directory.
# 2. Create 'POS' directory under top directory, move XDATCAR and split_xdat.py
#    there and run split_xdat.py to split XDATCAR to POSCAR_*.
# 3. Prepare INCAR, KPOINTS and POTCAR under top directory.
# 4. Modify run_wfn.sh to your need.
# 5. Run run_wfn.sh to generate the wave functions.

nmax=2000
ncore=20
vasp="vasp.mpi.5.4.1"

topdir=$(pwd)
for i in $(seq $nmax); do
        # create directory
        postfix=$(echo $i | awk '{printf "%04d", $1}')
        mkdir $topdir/$postfix && cd $topdir/$postfix

        # link files here
        ln -s ../INCAR .
        ln -s ../KPOINTS .
        ln -s ../POTCAR .
        ln -s ../POS/POSCAR_$postfix POSCAR
        
        # run vasp
        mpirun -np $ncore $vasp &> log
        wait
done
