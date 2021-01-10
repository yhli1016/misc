include(defs.m4)dnl
changecom("/*", "*/")dnl
#!/bin/bash
#BSUB -n NCPU 
#BSUB -W TIME:00
#BSUB -J NAME 
#BSUB -o log 

# 1) set the obname (BSUB -J)
# 2) set scratchroot
# 3) set number of run 
# 4) set ncores / time
# 5) adjust in INCAR  < 8 CPU --> NPAR=1, 16 CPU --> NPAR=2, >= 32 CPU --> NPAR=4  

name=${LSB_JOBNAME}
run=RUN

scratchroot=/cluster/scratch/zhangwenj
scratch=${scratchroot}/${name}

ifdef([RESTART], [
# for restarting
cp CONTCAR POSCAR])
cp POSCAR POSCAR_run${run}

mkdir -p ${scratch}

for i in INCAR POSCAR KPOINTS POTCAR vdw_kernel.bindat
do
cp ${i} ${scratch}
done

cd ${scratch}
mpirun vasp-54 > ${LS_SUBCWD}/out_run${run}

for i in OUTCAR CONTCAR
do
cp ${scratch}/${i} ${LS_SUBCWD}/
done

cp ${scratch}/vasprun.xml ${LS_SUBCWD}/vasprun_${run}.xml
