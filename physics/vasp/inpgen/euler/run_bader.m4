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
if [ "$istart" -ne 0 ]; then
    cp CONTCAR POSCAR
fi
cp POSCAR POSCAR_run${run}

for i in INCAR POSCAR KPOINTS POTCAR vdw_kernel.bindat
do
    cp ${i} ${scratch_dir}
done

# Run vasp
cd ${scratch_dir}
mpirun vasp_std > ${SLURM_SUBMIT_DIR}/out_run${run}

# Copy results back
for i in OUTCAR CONTCAR CHGCAR AECCAR0 AECCAR2 OSZICAR
do
    cp ${scratch_dir}/${i} ${SLURM_SUBMIT_DIR}/
done
cp ${scratch_dir}/vasprun.xml ${SLURM_SUBMIT_DIR}/vasprun_${run}.xml
