include(defs.m4)dnl
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

#------------------------------- Common header ---------------------------------
name=${LSB_JOBNAME}
run=RUN

scratchroot=/cluster/scratch/zhangwenj
scratch=${scratchroot}/${name}

mkdir -p ${scratch}

#--------------------------- Task-dependent scripts ----------------------------
# Copy files to scratch
## Shared input for all the images
for i in INCAR KPOINTS POTCAR vdw_kernel.bindat
do
cp ${i} ${scratch}
done

## Structure (POSCAR/CONTCAR) for each image
## For the 1st run, copy directories 00-0N to scratch.
## If restarting, only update the POSCARs.
for i in DIR_TOT
do
ifdef([RESTART], [dnl
# for restarting
cp ${i}/CONTCAR ${i}/POSCAR
cp ${i}/POSCAR ${i}/POSCAR_${run}
cp ${i}/POSCAR ${scratch}/${i}/POSCAR], [dnl
# for the first run
cp -r ${i} ${scratch}])
done

## OUTCAR of initial and final states for reference
for i in 00 NMAX
do
cp OUTCAR_$i $i/OUTCAR
cp OUTCAR_$i ${scratch}/$i/OUTCAR
done

# Run vasp
cd ${scratch}
mpirun vasp-54 > ${LS_SUBCWD}/out_run${run}

# Copy results back
for i in DIR_TS
do
cp ${scratch}/${i}/OUTCAR ${LS_SUBCWD}/${i}/OUTCAR
cp ${scratch}/${i}/CONTCAR ${LS_SUBCWD}/${i}/CONTCAR
done
cp ${scratch}/vasprun.xml ${LS_SUBCWD}/vasprun_${run}.xml
