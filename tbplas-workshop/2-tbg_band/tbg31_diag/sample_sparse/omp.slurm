#!/bin/bash
#SBATCH --partition=cpu_2d
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=tbg_tbpm
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Load prerequisites
module load compiler mkl mpi
module load anaconda3
export PYTHONPATH=$HOME/soft/tbplas

# Set number of threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Change to working directory and run the job
cd $SLURM_SUBMIT_DIR
mpirun -genv FI_PROVIDER=mlx python ./run_tbg.py
