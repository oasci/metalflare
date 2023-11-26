#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export PDB_ID="1JC0"
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=20
export METALFLARE_LOG_FILE_PATH="protein-prep.log"
METALFLARE_SAVE_DIR="structures/protein"

# Cleanup files from previous run
rm -rf $METALFLARE_SAVE_DIR
mkdir -p $METALFLARE_SAVE_DIR
rm -f $METALFLARE_LOG_FILE_PATH

# Get PDB file
wget https://files.rcsb.org/download/$PDB_ID.pdb -O $METALFLARE_SAVE_DIR/0-$PDB_ID.pdb

# Process PDB file
metalflare-select-atoms $METALFLARE_SAVE_DIR/0-$PDB_ID.pdb $METALFLARE_SAVE_DIR/1-$PDB_ID-chain-A.pdb --select_str chainID A and not resname HOH
metalflare-filter-pdb $METALFLARE_SAVE_DIR/1-$PDB_ID-chain-A.pdb  --output $METALFLARE_SAVE_DIR/2-$PDB_ID-filtered.pdb
metalflare-unify-resids $METALFLARE_SAVE_DIR/2-$PDB_ID-filtered.pdb --output $METALFLARE_SAVE_DIR/3-$PDB_ID-resid-fixes.pdb
metalflare-center $METALFLARE_SAVE_DIR/3-$PDB_ID-resid-fixes.pdb --output $METALFLARE_SAVE_DIR/4-$PDB_ID-centered.pdb
metalflare-rename-resname $METALFLARE_SAVE_DIR/4-$PDB_ID-centered.pdb MSE MET --output $METALFLARE_SAVE_DIR/5-$PDB_ID-resnames.pdb

pdb2pqr --log-level INFO --ff=AMBER --keep-chain --ffout=AMBER $METALFLARE_SAVE_DIR/5-$PDB_ID-resnames.pdb $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
cat $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.log >> $METALFLARE_LOG_FILE_PATH
rm $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.log
metalflare-rename-resname $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb HID HIS --output $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
metalflare-rename-resname $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb HIE HIS --output $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
metalflare-rename-resname $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb HIP HIS --output $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
metalflare-rename-resname $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb HOH WAT --output $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
metalflare-rename-resname $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb TIP WAT --output $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
metalflare-rename-resname $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb TIP3 WAT --output $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
metalflare-rename-resname $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb CYX CYS --output $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
metalflare-merge-pdbs $METALFLARE_SAVE_DIR/5-$PDB_ID-resnames.pdb $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb --output $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb

pdb4amber -i $METALFLARE_SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb > $METALFLARE_SAVE_DIR/7-$PDB_ID-pdb4amber.pdb 2> pdb4amber.err
cat pdb4amber.err >> $METALFLARE_LOG_FILE_PATH
rm pdb4amber.err

cp $METALFLARE_SAVE_DIR/7-$PDB_ID-pdb4amber.pdb $METALFLARE_SAVE_DIR/$PDB_ID-final.pdb

export METALFLARE_LOG=False
)
