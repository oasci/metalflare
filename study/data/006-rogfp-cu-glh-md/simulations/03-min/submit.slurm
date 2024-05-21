#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cluster=mpi
#SBATCH --partition=mpi
#SBATCH --account=jdurrant
#SBATCH --job-name=metalflare/003/03-min

# Load environment
module purge
module load gcc/8.2.0
module load openmpi/4.0.3
module load amber/22

cp ../02-prep/mol.inpcrd ./inputs
cp ../02-prep/mol.prmtop ./inputs

mkdir -p outputs/
./run.sh

touch CALC_DONE
