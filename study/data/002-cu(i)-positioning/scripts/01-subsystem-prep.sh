#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=10
export METALFLARE_LOG_FILE_PATH="01-subsystem-prep.log"

ORIGIN_DIR="../../001-rogfp-md/simulations/02-prep"

SAVE_DIR="../structures/01-initial"

# Cleanup files from previous run
rm -rf $SAVE_DIR
mkdir -p $SAVE_DIR
rm -f $METALFLARE_LOG_FILE_PATH

cp $INITIAL_PROTEIN $SAVE_DIR
./select-residues.py $ORIGIN_DIR/mol.prmtop $ORIGIN_DIR/mol.inpcrd $SAVE_DIR/subsystem.xyz --resi 145 202 --within 6.0

export METALFLARE_LOG=False
)
