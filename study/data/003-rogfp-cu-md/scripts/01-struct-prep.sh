#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=10
export METALFLARE_LOG_FILE_PATH="01-struct-prep.log"

rm $METALFLARE_LOG_FILE_PATH
rm $PDB_PATH

PDB_PATH="../structures/protein/1JC0-Cu.pdb"

python3 ./01-integrate-copper.py
metalflare-unify-resids $PDB_PATH --output $PDB_PATH

export METALFLARE_LOG=False
)
