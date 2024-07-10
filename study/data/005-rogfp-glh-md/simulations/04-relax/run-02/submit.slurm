#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=48
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=3-00:00:00
#SBATCH --cluster=mpi
#SBATCH --partition=mpi
#SBATCH --account=jdurrant
#SBATCH --job-name=metalflare/005/04-relax/run-02

# Load environment
module purge
module load gcc/8.2.0
module load openmpi/4.0.3
module load amber/22

cp ../../02-prep/mol.inpcrd ./inputs
cp ../../02-prep/mol.prmtop ./inputs

mkdir -p outputs/
./run.sh

touch CALC_DONE
