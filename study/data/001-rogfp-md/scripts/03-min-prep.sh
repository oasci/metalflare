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
COORD_PATH="$SIMULATIONS_DIR/02-prep/mol.inpcrd"
TOPO_PATH="$SIMULATIONS_DIR/02-prep/mol.prmtop"

SAVE_DIR="$SIMULATIONS_DIR/03-min"
INPUT_DIR="$SAVE_DIR/inputs"
OUTPUT_DIR="$SAVE_DIR/outputs"
RUN_PATH="$SAVE_DIR/run.sh"
SLURM_PATH="$SAVE_DIR/submit.slurm"
JOB_NAME="metalflare/001/03-min"


# Cleanup files from previous run
rm -rf $SAVE_DIR
mkdir -p $SAVE_DIR
rm -rf $INPUT_DIR
mkdir -p $INPUT_DIR
rm -rf $OUTPUT_DIR
mkdir -p $OUTPUT_DIR
rm -f $METALFLARE_LOG_FILE_PATH

./prep_sim.py $TOPO_PATH $COORD_PATH $INPUT_DIR $RUN_PATH $JOB_NAME $SLURM_PATH \
--yaml $SIMULATIONS_DIR/base.yml $SIMULATIONS_DIR/03-min.yml $SIMULATIONS_DIR/slurm.yml

cp $COORD_PATH $INPUT_DIR/mol.inpcrd
cp $TOPO_PATH $INPUT_DIR/mol.prmtop

export METALFLARE_LOG=False
)
