#!/bin/bash
#SBATCH --account=def-sperling
#SBATCH --mail-user=ccid@ualberta.ca
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# sample slurm commands

# Run script type
# sbatch cedar_sbatch.sh

# List queued jobs
# squeue -u $USER

# Cancel all jobs
# scancel -u $USER

# Add commands below

