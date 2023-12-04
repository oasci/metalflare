#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=10
export METALFLARE_LOG_FILE_PATH="05-prod-prep.log"

SIMULATIONS_DIR="../simulations"
COORD_PATH="$SIMULATIONS_DIR/02-prep/mol.inpcrd"
TOPO_PATH="$SIMULATIONS_DIR/02-prep/mol.prmtop"

REPLICATES=3

for ((i=1; i<=$REPLICATES; i++)); do
    suffix=$(printf "%02d" "$i")

    export RUN_NAME="run-$suffix"
    SAVE_DIR="$SIMULATIONS_DIR/05-prod/$RUN_NAME"
    INPUT_DIR="$SAVE_DIR/inputs"
    OUTPUT_DIR="$SAVE_DIR/outputs"
    RUN_PATH="$SAVE_DIR/run.sh"
    SLURM_PATH="$SAVE_DIR/submit.slurm"
    JOB_NAME="metalflare/001/05-prod/$RUN_NAME"


    # Cleanup files from previous run
    rm -rf $SAVE_DIR
    mkdir -p $SAVE_DIR
    rm -rf $INPUT_DIR
    mkdir -p $INPUT_DIR
    rm -rf $OUTPUT_DIR
    mkdir -p $OUTPUT_DIR
    rm -f $METALFLARE_LOG_FILE_PATH

    temp_yaml="$(mktemp)"
    cp $SIMULATIONS_DIR/05-prod.yml $temp_yaml
    sed -i "s/{{ RUN_NAME }}/$RUN_NAME/g" $temp_yaml

    ./prep_sim.py $TOPO_PATH $COORD_PATH $INPUT_DIR $RUN_PATH $JOB_NAME $SLURM_PATH \
    --yaml $SIMULATIONS_DIR/base.yml $SIMULATIONS_DIR/slurm.yml $temp_yaml

    cp $COORD_PATH $INPUT_DIR/mol.inpcrd
    cp $TOPO_PATH $INPUT_DIR/mol.prmtop
done

export METALFLARE_LOG=False
)
