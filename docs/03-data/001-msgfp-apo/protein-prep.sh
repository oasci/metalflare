#!/usr/bin/env bash

(
cd "$(dirname "$0")"

export PDB_ID="1JC0"
export METALFLARE_LOG=True
export METALFLARE_STDOUT=False
export METALFLARE_LOG_LEVEL=20
export METALFLARE_LOG_FILE_PATH="prep.log"

rm -f $PDB_ID.pdb
wget https://files.rcsb.org/download/$PDB_ID.pdb -O structures/0-$PDB_ID.pdb

mkdir -p structures
rm -f prep.log

metalflare-filter-pdb structures/0-$PDB_ID.pdb --output structures/1-$PDB_ID-filtered.pdb
metalflare-unify-resids structures/1-$PDB_ID-filtered.pdb --output structures/2-$PDB_ID-resids.pdb

export METALFLARE_LOG=False
)
