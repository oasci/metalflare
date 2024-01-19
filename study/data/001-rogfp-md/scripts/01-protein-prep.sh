#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export PDB_ID="1JC0"
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=10
export METALFLARE_LOG_FILE_PATH="01-protein-prep.log"
SAVE_DIR="../structures/protein"

# Cleanup files from previous run
rm -rf $SAVE_DIR
mkdir -p $SAVE_DIR
rm -f $METALFLARE_LOG_FILE_PATH

# Get PDB file
wget https://files.rcsb.org/download/$PDB_ID.pdb -O $SAVE_DIR/0-$PDB_ID.pdb

# Process PDB file
metalflare-select-atoms $SAVE_DIR/0-$PDB_ID.pdb $SAVE_DIR/1-$PDB_ID-chain-A.pdb --select_str "chainID A"
metalflare-filter-pdb $SAVE_DIR/1-$PDB_ID-chain-A.pdb  --output $SAVE_DIR/2-$PDB_ID-filtered.pdb
metalflare-center $SAVE_DIR/2-$PDB_ID-filtered.pdb --output $SAVE_DIR/3-$PDB_ID-centered.pdb
metalflare-minimize-box $SAVE_DIR/3-$PDB_ID-centered.pdb --output $SAVE_DIR/4-$PDB_ID-rotated.pdb

metalflare-unify-resids $SAVE_DIR/4-$PDB_ID-rotated.pdb --output $SAVE_DIR/5-$PDB_ID-residues.pdb
metalflare-rename-resname $SAVE_DIR/5-$PDB_ID-residues.pdb MSE MET --output $SAVE_DIR/5-$PDB_ID-residues.pdb
metalflare-rename-resname $SAVE_DIR/5-$PDB_ID-residues.pdb CYS CYM --include 145 202 --output $SAVE_DIR/5-$PDB_ID-residues.pdb

pdb2pqr --log-level INFO --ff=AMBER --keep-chain --ffout=AMBER $SAVE_DIR/5-$PDB_ID-residues.pdb $SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
cat $SAVE_DIR/6-$PDB_ID-pdb2pqr.log >> $METALFLARE_LOG_FILE_PATH
rm $SAVE_DIR/6-$PDB_ID-pdb2pqr.log
metalflare-merge-pdbs $SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb $SAVE_DIR/5-$PDB_ID-residues.pdb --output $SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb

# Put all residue renames here.
metalflare-rename-resname $SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb HOH WAT --output $SAVE_DIR/7-$PDB_ID-resnames.pdb
metalflare-rename-resname $SAVE_DIR/7-$PDB_ID-resnames.pdb TIP WAT --output $SAVE_DIR/7-$PDB_ID-resnames.pdb
metalflare-rename-resname $SAVE_DIR/7-$PDB_ID-resnames.pdb TIP3 WAT --output $SAVE_DIR/7-$PDB_ID-resnames.pdb

pdb4amber -i $SAVE_DIR/7-$PDB_ID-resnames.pdb > $SAVE_DIR/8-$PDB_ID-pdb4amber.pdb 2> pdb4amber.err
cat pdb4amber.err >> $METALFLARE_LOG_FILE_PATH
rm pdb4amber.err

cp $SAVE_DIR/8-$PDB_ID-pdb4amber.pdb $SAVE_DIR/$PDB_ID-final.pdb

export METALFLARE_LOG=False
)
