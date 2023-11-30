#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=3-00:00:00
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
#SBATCH --account=jdurrant

# Load environment
module purge
module load gcc/8.2.0
module load openmpi/4.0.3
module load amber/22

total_n_cores=$(($SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE))
