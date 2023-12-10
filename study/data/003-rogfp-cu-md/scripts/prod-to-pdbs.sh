#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

SIMULATIONS_DIR="../simulations"
REPLICATES=1

for ((i=1; i<=$REPLICATES; i++)); do
    suffix=$(printf "%02d" "$i")

    RUN_NAME="run-$suffix"

    STAGE_DIR="$SIMULATIONS_DIR/05-prod/$RUN_NAME"
    OUTPUT_DIR="$STAGE_DIR/outputs"
    COORD_PATH="$OUTPUT_DIR/08_prod_npt.nc"
    TOPO_PATH="$SIMULATIONS_DIR/02-prep/mol.prmtop"

    metalflare-pdb $OUTPUT_DIR/08_prod_npt_raw.pdb --files $TOPO_PATH $COORD_PATH --stride 250 --select not resname WAT and not resname Cl- and not resname Na+
    metalflare-pdb-align $OUTPUT_DIR/08_prod_npt_raw.pdb $OUTPUT_DIR/08_prod_npt.pdb --select protein
    rm $OUTPUT_DIR/08_prod_npt_raw.pdb
done

)
