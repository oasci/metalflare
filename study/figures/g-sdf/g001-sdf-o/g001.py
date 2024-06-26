from collections.abc import Mapping

import numpy as np
import pymol
from pymol import cmd, util

# https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha

ISO_VALUE = 0.0014


def move_object(mobile_object: str, x: Mapping[float]) -> None:
    """Move PyMOL object's display matrix. Transformations are made in PyMOLs internal
    space and are not physical changes.

    Args:
        mobile_object: PyMOL object name to be transformed.
        x: Custom transformation vector of
            `[tran_x, trans_y, trans_z, rot_x, rot_y, rot_z]`.
    """
    translate = x[:3]
    angles = x[3:]
    if isinstance(translate, np.ndarray):
        translate = translate.tolist()
    cmd.translate(translate, object=mobile_object, camera=0)
    cmd.rotate("x", angles[0], object=mobile_object, camera=0)
    cmd.rotate("y", angles[1], object=mobile_object, camera=0)
    cmd.rotate("z", angles[2], object=mobile_object, camera=0)


# Start PyMOL session
pymol.finish_launching()

# Setup colors
cmd.bg_color("white")
cmd.set_color("red_color", [30, 46, 121])
cmd.set_color("oxd_color", [236, 64, 103])
cmd.set_color("cu_color", [249, 151, 82])

# Load files
cmd.load(
    "../../../analysis/001-rogfp-md/data/structure/average-structure.pdb", "red_protein"
)
cmd.load(
    "../../../analysis/004-rogfp-oxd-md/data/structure/average-structure.pdb",
    "oxd_protein",
)
cmd.load(
    "../../../analysis/003-rogfp-cu-md/data/structure/average-structure.pdb",
    "cu_protein",
)

# Display specific residues as sticks
cmd.select("proton_wire", "(resi 201 or resi 203 or resi 220)")
cmd.show("sticks", "proton_wire")
cmd.select("close_residues", "(resi 143 or resi 146)")
cmd.show("sticks", "close_residues")
cmd.set("sphere_transparency", 0.9, "close_residues")
cmd.set("stick_transparency", 0.9, "close_residues")
cmd.select("cys", "(resi 145 or resi 202)")
cmd.hide("cartoon", "(not proton_wire and not cys)")
cmd.select("cro", "(resn CRO)")
cmd.set("sphere_transparency", 0.3, "cro")
cmd.set("stick_transparency", 0.3, "cro")

# Load density maps
cmd.load("../../../analysis/001-rogfp-md/data/sdf/resid65_oh-o.dx", "red")
cmd.load("../../../analysis/004-rogfp-oxd-md/data/sdf/resid65_oh-o.dx", "oxd")
cmd.load("../../../analysis/003-rogfp-cu-md/data/sdf/resid65_oh-o.dx", "cu")

# Generate and adjust iso surfaces
cmd.isosurface("red_iso", "red", 1.0)
cmd.isosurface("oxd_iso", "oxd", 1.0)
cmd.isosurface("cu_iso", "cu", 1.0)

cmd.isolevel("red_iso", ISO_VALUE)
cmd.isolevel("oxd_iso", ISO_VALUE)
cmd.isolevel("cu_iso", ISO_VALUE)

cmd.color("red_color", "red_iso")
cmd.color("oxd_color", "oxd_iso")
cmd.color("cu_color", "cu_iso")

cmd.set("transparency", 0.2, "red_iso")
cmd.set("transparency", 0.4, "oxd_iso")
cmd.set("transparency", 0.3, "cu_iso")

# Load and perform alignment transformations
oxd_transform = np.load("transform_oxd.npy").tolist()
cu_transform = np.load("transform_cu.npy").tolist()

move_object("oxd_protein", oxd_transform)
move_object("cu_protein", cu_transform)
move_object("oxd_iso", oxd_transform)
move_object("cu_iso", cu_transform)

# Disable other objects
cmd.disable("oxd_protein")
cmd.disable("cu_protein")
# cmd.disable("oxd_iso")

# Labels
# cmd.label(selection="model red_protein and resi 220 and name CD", expression="resn")
# cmd.label(selection="model red_protein and resi 201 and name OG1", expression="resn")
# cmd.label(selection="model red_protein and resi 203 and name OG", expression="resn")
cmd.set("label_font_id", 7)
cmd.set("label_size", 22)

# Set view
cmd.color("grey70", "element C")
cmd.set_view(
    """
(\
     0.115785033,   -0.774704337,   -0.621626854,\
     0.524430394,    0.579164386,   -0.624106050,\
     0.843521476,   -0.253746152,    0.473351479,\
    -0.002147511,   -0.000185485,  -27.974308014,\
    38.576263428,   35.590415955,   39.145816803,\
  -150.627853394,  206.387130737,  -20.000000000 )
"""
)

cmd.refresh()

# Figure settings
cmd.set("ray_trace_fog", 0)
cmd.set("depth_cue", 0)
util.performance(0)
cmd.set("antialias", 4)
cmd.set("hash_max", 300)
cmd.set("ray_trace_mode", 1)

cmd.enable("oxd_iso")
cmd.disable("cu_iso")
cmd.png("red_oxd-sdf-o.png", dpi=1000)
cmd.disable("oxd_iso")
cmd.enable("cu_iso")
cmd.png("red_cu-sdf-o.png", dpi=1000)
cmd.refresh()
