#!/bin/bash
#BSUB -n <NCPU>
#BSUB -W <TIME>:00
#BSUB -J <NAME>
#BSUB -o log

# 1) set the obname (BSUB -J)
# 2) set scratchroot
# 3) set number of run
# 4) set ncores / time
# 5) adjust in INCAR  < 8 CPU --> NPAR=1, 16 CPU --> NPAR=2, >= 32 CPU --> NPAR=4

#------------------------------- Common header ---------------------------------
name=${LSB_JOBNAME}
run=<RUN>
istart=<ISTART>

scratchroot=/cluster/scratch/zhangwenj
scratch=${scratchroot}/${name}

mkdir -p ${scratch}

#--------------------------- Task-dependent scripts ----------------------------
# Copy files to scratch
if [ "$istart" -ne 0 ]; then
    cp CONTCAR POSCAR
fi
cp POSCAR POSCAR_run${run}

for i in INCAR POSCAR KPOINTS POTCAR vdw_kernel.bindat
do
cp ${i} ${scratch}
done

# Run vasp
cd ${scratch}
mpirun vasp-54 > ${LS_SUBCWD}/out_run${run}

# Copy results back
for i in OUTCAR CONTCAR CHGCAR AECCAR0 AECCAR2
do
cp ${scratch}/${i} ${LS_SUBCWD}/
done
cp ${scratch}/vasprun.xml ${LS_SUBCWD}/vasprun_${run}.xml
