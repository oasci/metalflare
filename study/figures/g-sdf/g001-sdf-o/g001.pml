
# Setup colors
set_color red_color, [30, 46, 121]
set_color oxd_color, [236, 64, 103]
set_color cu_color, [249, 151, 82]

# Load files
load ../../../analysis/001-rogfp-md/data/sdf/resid65_oh-o.dx, red
load ../../../analysis/004-rogfp-oxd-md/data/sdf/resid65_oh-o.dx, oxd
load ../../../analysis/003-rogfp-cu-md/data/sdf/resid65_oh-o.dx, cu

# Generate surfaces
isomesh red_iso, red
isomesh oxd_iso, oxd
isomesh cu_iso, cu

# Fix surface levels, 0.05
isolevel red_iso, 0.1
isolevel oxd_iso, 0.1
isolevel cu_iso, 0.1

# Change colors
color red_color, red_iso
color oxd_color, oxd_iso
color cu_color, cu_iso

# Update origin
center red_iso

# Translate
# oxd -> red  [0.648, 4.086, -1.367]
# cu -> red   [-1.07, 4.712, 2.017]
pseudoatom red_cro_oh, pos=[42.554, 38.993, 31.583]
pseudoatom oxd_cro_oh, pos=[41.906, 34.907, 32.950]
pseudoatom cu_cro_oh, pos=[43.624, 34.281, 29.566]
