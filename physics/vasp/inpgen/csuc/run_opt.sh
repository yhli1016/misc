#! /bin/bash
#SBATCH --ntasks=<NCPU>
#SBATCH --cpus-per-task=1
#SBATCH --time=<TIME>:00:00
#SBATCH --job-name=<NAME>
#SBATCH --partition=std
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

scratchroot=/scratch/wzhang
scratch=${scratchroot}/${SLURM_JOB_NAME}

module load vasp/5.4.4_beef_vtst
export PATH=/prod/apps/vasp/std/5.4.4_vtst_beef/bin:$PATH
ulimit -s unlimited

# Remove existing scratch directory if istart == 0
# We must check $scratchroot carefully to avoid removing important files and
# directories by accident.
if [ "$scratchroot" == "$scratch" ]; then
    echo "ERROR: empty job name"
    exit 1
else
    if [ "$istart" -eq 0 ]; then
        rm -r ${scratch}
    fi
fi
mkdir -p ${scratch}

############################# Task-dependent part ##############################
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
srun vasp_std > ${SLURM_SUBMIT_DIR}/out_run${run}

# Copy results back
for i in OUTCAR CONTCAR OSZICAR
do
    cp ${scratch}/${i} ${SLURM_SUBMIT_DIR}/
done
cp ${scratch}/vasprun.xml ${SLURM_SUBMIT_DIR}/vasprun_${run}.xml
