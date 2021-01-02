include(defs.m4)dnl
changecom("/*", "*/")dnl
#!/bin/bash
#BSUB -n NCPU
#BSUB -W TIME:00
#BSUB -J NAME
#BSUB -o log 

# 1) set the jobname (BSUB -J)
# 2) set scratchroot
# 3) set number of run 
# 4) set ncores / time
# 5) adjust in INCAR  < 8 CPU --> NPAR=1, 16 CPU --> NPAR=2, >= 32 CPU --> NPAR=4  

name=${LSB_JOBNAME}
run=RUN

scratchroot=/cluster/scratch/wenjzhang
scratch=${scratchroot}/${name}

mkdir -p ${scratch}

for i in INCAR KPOINTS POTCAR vdw_kernel.bindat
do
cp ${i} ${scratch}
done

ifdef([RESTART], [
# for restarting
for i in DIR_TOT
do
cp ${i}/CONTCAR ${i}/POSCAR
cp ${i}/POSCAR ${i}/POSCAR_${run}
cp ${i}/POSCAR ${scratch}/${i}/POSCAR
done], [
# for the first run
for i in DIR_TOT
do
cp -r ${i} ${scratch}
done])

for i in 00 NMAX
do
cp OUTCAR_$i ${scratch}/$i/OUTCAR
cp OUTCAR_$i $i/OUTCAR
done

cd ${scratch}
mpirun vasp-54 > ${LS_SUBCWD}/out_run${run}

for i in DIR_TS
do
cp ${scratch}/${i}/OUTCAR ${LS_SUBCWD}/${i}/OUTCAR
cp ${scratch}/${i}/CONTCAR ${LS_SUBCWD}/${i}/CONTCAR
done

cp ${scratch}/vasprun.xml ${LS_SUBCWD}/vasprun_${run}.xml
