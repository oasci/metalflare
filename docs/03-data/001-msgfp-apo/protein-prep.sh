#!/usr/bin/env bash
(
cd "$(dirname "$0")"

export PDB_ID="1JC0"

rm -f $PDB_ID.pdb
wget https://files.rcsb.org/download/$PDB_ID.pdb

metalflare-filter-pdb $PDB_ID.pdb --output $PDB_ID-filtered.pdb

)
