#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=10
export METALFLARE_LOG_FILE_PATH="03-min-prep.log"

SIMULATIONS_DIR="../simulations"
COORD_PATH="$SIMULATIONS_DIR/02-amber-prep/mol.inpcrd"
TOPO_PATH="$SIMULATIONS_DIR/02-amber-prep/mol.prmtop"
SLURM_PATH="$SIMULATIONS_DIR/02-amber-prep/run.slurm"
SAVE_DIR="$SIMULATIONS_DIR/03-min"


# Cleanup files from previous run
rm -rf $SAVE_DIR
mkdir -p $SAVE_DIR
rm -f $METALFLARE_LOG_FILE_PATH

./prep_sim.py $TOPO_PATH $COORD_PATH $SLURM_PATH --yaml $SIMULATIONS_DIR/base.yml $SIMULATIONS_DIR/slurm.yml

export METALFLARE_LOG=False
)
