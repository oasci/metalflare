#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=10
export METALFLARE_LOG_FILE_PATH="02-amber-prep.log"

PDB_PATH="../structures/protein/1JC0-final.pdb"
YAML_PATH="../simulations/base.yml"

SAVE_DIR="../simulations/02-prep"
TOPO_PATH="$SAVE_DIR/mol.prmtop"
COORD_PATH="$SAVE_DIR/mol.inpcrd"



# Cleanup files from previous run
rm -rf "$SAVE_DIR"
mkdir "$SAVE_DIR"
rm -f $METALFLARE_LOG_FILE_PATH

metalflare-validate-context $YAML_PATH metalflare.simulation.amber.contexts.AmberContextValidator
./run_tleap.py $PDB_PATH $TOPO_PATH $COORD_PATH --yaml $YAML_PATH

export METALFLARE_LOG=False
)
