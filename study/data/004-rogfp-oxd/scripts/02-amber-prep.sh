#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=10
export METALFLARE_LOG_FILE_PATH="02-amber-prep.log"
METHODS_PATH="../../../methods"

PDB_PATH="../structures/protein/1JC1-final.pdb"
YAML_PATH="$METHODS_PATH/04-amber-simulations/base.yml"

SAVE_DIR="../simulations/02-prep"
TOPO_PATH="$SAVE_DIR/mol.prmtop"
COORD_PATH="$SAVE_DIR/mol.inpcrd"
OUTPUT_PDB_PATH="$SAVE_DIR/mol.pdb"



# Cleanup files from previous run
rm -rf "$SAVE_DIR"
mkdir -p "$SAVE_DIR"
rm -f $METALFLARE_LOG_FILE_PATH

metalflare-validate-context $YAML_PATH metalflare.simulation.amber.contexts.AmberContextValidator
metalflare-tleap $PDB_PATH $TOPO_PATH $COORD_PATH --yaml $YAML_PATH --work_dir "$(dirname "$0")"
metalflare-pdb $OUTPUT_PDB_PATH --files $TOPO_PATH $COORD_PATH

export METALFLARE_LOG=False
)
