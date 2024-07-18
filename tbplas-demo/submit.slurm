#!/bin/bash
#SBATCH --partition=cpu_2d
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=si_dc
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Load prerequisites
conda deactivate
module load compiler mkl mpi
module load anaconda3
export PYTHONPATH=$HOME/soft/tbplas-demo:$HOME/soft/tbplas-develop

# Set number of threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Change to working directory and run the job
cd $SLURM_SUBMIT_DIR
mpirun -genv FI_PROVIDER=mlx python $HOME/soft/tbplas-demo/main.py
