#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=10
export METALFLARE_LOG_FILE_PATH="amber-prep.log"
METALFLARE_SAVE_DIR="simulations/prep"


# Cleanup files from previous run
rm -rf $METALFLARE_SAVE_DIR
mkdir -p $METALFLARE_SAVE_DIR
rm -f $METALFLARE_LOG_FILE_PATH

metalflare-validate-context simulations/base.yml metalflare.simulation.amber.contexts.AmberContextValidator
./run_tleap.py structures/protein/1JC0-final.pdb --yaml simulations/base.yml

export METALFLARE_LOG=False
)
