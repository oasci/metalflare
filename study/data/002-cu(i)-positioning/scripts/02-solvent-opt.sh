#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

SAVE_DIR="../calculations/01-solvent-opt"
INPUT_FILE="$SAVE_DIR/opt.in"

# Cleanup files from previous run
rm -rf $SAVE_DIR
mkdir -p $SAVE_DIR

cat > $INPUT_FILE << EOF
\$constrain
   force constant=0.5
   elements: C,N,O,S
\$end
EOF

cd $SAVE_DIR
INITIAL_XYZ="../../structures/01-initial/subsystem.xyz"
xtb $INITIAL_XYZ --opt loose --gfn 2 --alpb water --input opt.in > opt.log
)
