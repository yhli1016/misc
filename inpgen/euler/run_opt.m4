include(common.m4)dnl

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
