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
METHODS_DIR="../../../methods"
YAML_DIR="$METHODS_DIR/04-amber-simulations"

COORD_PATH="$SIMULATIONS_DIR/02-prep/mol.inpcrd"
TOPO_PATH="$SIMULATIONS_DIR/02-prep/mol.prmtop"

REPLICATES=3

for ((i=1; i<=$REPLICATES; i++)); do
    suffix=$(printf "%02d" "$i")

    RUN_NAME="run-$suffix"

    STAGE_DIR="$SIMULATIONS_DIR/05-prod/$RUN_NAME"
    INPUT_DIR="$STAGE_DIR/inputs"
    OUTPUT_DIR="$STAGE_DIR/outputs"

    JOB_NAME="metalflare/001/05-prod/$RUN_NAME"
    WRITE_DIR="$INPUT_DIR"
    RUN_PATH="$STAGE_DIR/run.sh"
    SLURM_PATH="$STAGE_DIR/submit.slurm"
    PREP_CLASS_STRING="metalflare.simulation.amber.run.AmberRunPrep"

    # Cleanup files from previous run
    rm -rf $SAVE_DIR
    mkdir -p $SAVE_DIR
    rm -rf $INPUT_DIR
    mkdir -p $INPUT_DIR
    rm -rf $OUTPUT_DIR
    mkdir -p $OUTPUT_DIR
    rm -f $METALFLARE_LOG_FILE_PATH

    temp_yaml="$(mktemp)"
    cp $YAML_DIR/05-prod.yml $temp_yaml
    sed -i "s/{{ RUN_NAME }}/$RUN_NAME/g" $temp_yaml

    metalflare-prep-sims $JOB_NAME $WRITE_DIR $RUN_PATH $SLURM_PATH $PREP_CLASS_STRING \
    --yaml $temp_yaml $YAML_DIR/slurm-runs.yml $YAML_DIR/base.yml

    chmod +x $RUN_PATH

done

export METALFLARE_LOG=False
)
