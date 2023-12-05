#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

SAVE_DIR="../calculations/02-dock-copper"
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

INITIAL_XYZ="../01-solvent-opt/xtbopt.xyz"
cp $INITIAL_XYZ initial.xyz
sed -i -e '$aCu          50.94271                33.39869                35.24655' initial.xyz
sed -i 's/314/315/' initial.xyz

CHRG_PATH="../../structures/01-initial/.CHRG"
cp $CHRG_PATH .CHRG
current_number=$(<.CHRG)
new_number=$((current_number + 1))
echo "$new_number" > .CHRG

xtb initial.xyz --opt loose --gfn 2 --alpb water --input opt.in > opt.log
)
