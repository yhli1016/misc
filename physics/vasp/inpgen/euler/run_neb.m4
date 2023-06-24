#! /bin/bash
#SBATCH --ntasks=<NCPU>
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=<MEM>
#SBATCH --time=<TIME>:00:00
#SBATCH --job-name=<NAME>
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# NOTES:
# [1] set the number of cores (--ntasks)
# [2] set the wall time (--time)
# [3] set the job name (--job-name)
# [4] set the number of run and istart
# [5] adjust INCAR:
#     < 8 CPU: NPAR=1, 16 CPU: NPAR=2, >= 32 CPU: NPAR=4

################################# Common part ##################################
run=<RUN>
istart=<ISTART>

scratch_dir=${SCRATCH}/${SLURM_JOB_NAME}

# Remove existing scratch_dir directory if istart == 0
if [ "$SCRATCH" == "$scratch_dir" ]; then
    echo "ERROR: empty job name"
    exit 1
else
    if [ "$istart" -eq 0 ]; then
        rm -r ${scratch_dir}
    fi
fi
mkdir -p ${scratch_dir}

############################# Task-dependent part ##############################
# Copy files to scratch_dir
## Shared input for all the images
for i in INCAR KPOINTS POTCAR vdw_kernel.bindat
do
    cp ${i} ${scratch_dir}
done

## Structure (POSCAR/CONTCAR) for each image
## For the 1st run, copy directories 00-0N to scratch_dir.
## If restarting, only update the POSCARs.
for i in <DIR_TOT>
do
    if [ "$istart" -eq 0 ]; then
        cp -r ${i} ${scratch_dir}
    else
        cp ${i}/CONTCAR ${i}/POSCAR
        cp ${i}/POSCAR ${i}/POSCAR_${run}
        cp ${i}/POSCAR ${scratch_dir}/${i}/POSCAR
    fi
done

## OUTCAR of initial and final states for reference
for i in 00 <NMAX>
do
    cp OUTCAR_$i $i/OUTCAR
    cp OUTCAR_$i ${scratch_dir}/$i/OUTCAR
done

# Run vasp
cd ${scratch_dir}
mpirun vasp_std > ${SLURM_SUBMIT_DIR}/out_run${run}

# Copy results back
for i in <DIR_TS>
do
    cp ${scratch_dir}/${i}/OUTCAR ${SLURM_SUBMIT_DIR}/${i}/OUTCAR
    cp ${scratch_dir}/${i}/CONTCAR ${SLURM_SUBMIT_DIR}/${i}/CONTCAR
    cp ${scratch_dir}/${i}/OSZICAR ${SLURM_SUBMIT_DIR}/${i}/OSZICAR
done
cp ${scratch_dir}/vasprun.xml ${SLURM_SUBMIT_DIR}/vasprun_${run}.xml
