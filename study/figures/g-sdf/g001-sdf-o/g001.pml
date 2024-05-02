
# Setup colors
set_color red_color, [30, 46, 121]
set_color oxd_color, [236, 64, 103]
set_color cu_color, [249, 151, 82]

# Load files
load ../../../data/001-rogfp-md/structures/protein/1JC0-final.pdb, 1JC0
load ../../../analysis/001-rogfp-md/data/sdf/resid65_oh-o.dx, red
load ../../../analysis/004-rogfp-oxd-md/data/sdf/resid65_oh-o.dx, oxd
load ../../../analysis/003-rogfp-cu-md/data/sdf/resid65_oh-o.dx, cu

# Show atoms
select proton_wire, (resi 143 or resi 146 or resi 201 or resi 220)
show sticks, proton_wire
deselect
select cys, (resi 145 or resi 202)
show sticks, cys
deselect

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

# Translate densities
# oxd -> red  [0.648, 4.086, -1.367]
# cu -> red   [-1.07, 4.712, 2.017]
pseudoatom red_cro_oh, pos=[42.554, 38.993, 31.583]
# pseudoatom oxd_cro_oh, pos=[41.906, 34.907, 32.950]
# pseudoatom cu_cro_oh, pos=[43.624, 34.281, 29.566]

translate [0.648, 4.086, -1.367], object = oxd_iso, camera = 0
translate [-1.07, 4.712, 2.017], object = cu_iso, camera = 0

# Translate protein
# gfp -> red [7.468, 3.095, -0.106]
translate [35.086, 35.898, 31.689], object = 1JC0, camera = 0
#rotate x, -100, origin=[42.554, 38.993, 31.583], camera=0
#rotate y, -50, origin=[42.554, 38.993, 31.583], camera=0

select cro_oh, (resn cro and name OH)
deselect
center cro_oh
