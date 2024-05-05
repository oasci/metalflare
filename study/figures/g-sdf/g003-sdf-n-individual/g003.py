from collections.abc import Mapping

import numpy as np
import pymol
from pymol import cmd, util

# https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha

ISO_VALUE = 0.09


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
cmd.select("cys-group", "(resi 145 or resi 202 or resname CYM)")
cmd.hide("cartoon", "(not proton_wire and not cys-group)")
cmd.hide("everything", "(cys-group)")
cmd.select("cro", "(resn CRO)")
cmd.show("sticks", "cro")


map_data = {
    "red": {
        "cro": "../../../analysis/001-rogfp-md/data/sdf/resnamecro_n.dx",
        "tyr143": "../../../analysis/001-rogfp-md/data/sdf/resid143_n.dx",
        "his146": "../../../analysis/001-rogfp-md/data/sdf/resid146_n.dx",
        "thr201": "../../../analysis/001-rogfp-md/data/sdf/resid201_n.dx",
        "ser203": "../../../analysis/001-rogfp-md/data/sdf/resid203_n.dx",
        "glu220": "../../../analysis/001-rogfp-md/data/sdf/resid220_n.dx",
    },
    "oxd": {
        "cro": "../../../analysis/004-rogfp-oxd-md/data/sdf/resnamecro_n.dx",
        "tyr143": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid143_n.dx",
        "his146": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid146_n.dx",
        "thr201": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid201_n.dx",
        "ser203": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid203_n.dx",
        "glu220": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid220_n.dx",
    },
    "cu": {
        "cro": "../../../analysis/003-rogfp-cu-md/data/sdf/resnamecro_n.dx",
        "tyr143": "../../../analysis/003-rogfp-cu-md/data/sdf/resid143_n.dx",
        "his146": "../../../analysis/003-rogfp-cu-md/data/sdf/resid146_n.dx",
        "thr201": "../../../analysis/003-rogfp-cu-md/data/sdf/resid201_n.dx",
        "ser203": "../../../analysis/003-rogfp-cu-md/data/sdf/resid203_n.dx",
        "glu220": "../../../analysis/003-rogfp-cu-md/data/sdf/resid220_n.dx",
    },
}


# Load alignment transformations
trans_data = {
    "oxd": np.load("../g001-sdf-o/transform_oxd.npy").tolist(),
    "cu": np.load("../g001-sdf-o/transform_cu.npy").tolist(),
}
move_object("oxd_protein", trans_data["oxd"])
move_object("cu_protein", trans_data["cu"])

# Load density maps
for sys_label, map_info in map_data.items():
    for res_label, map_path in map_info.items():
        map_label = f"{sys_label}_{res_label}"
        cmd.load(map_path, map_label)
        if sys_label in trans_data.keys():
            move_object(map_label, trans_data[sys_label])
        cmd.group(name="maps", members=map_label, action="add")

        iso_label = map_label + "_iso"
        cmd.isosurface(name=iso_label, map=map_label)
        cmd.isolevel(iso_label, ISO_VALUE)
        cmd.group(name=f"{res_label}_iso", members=iso_label, action="add")

        cmd.color(f"{sys_label}_color", iso_label)
        cmd.set("transparency", 0.3, iso_label)


# Disable other objects
cmd.disable("oxd_protein")
cmd.disable("cu_protein")
cmd.disable("maps")

# Set view
cmd.color("grey70", "element C")

# Figure settings
util.performance(0)
cmd.space("cmyk")
cmd.set("ray_trace_fog", 0)
cmd.set("depth_cue", 0)
cmd.set("antialias", 4)
cmd.set("hash_max", 300)
cmd.set("ray_trace_mode", 1)
cmd.set("opaque_background", 0)
cmd.set("cartoon_discrete_colors", 0)
cmd.set("ambient", 0)

cmd.refresh()


TRANSPARENCY_VALUE = 0.8
cmd.enable("cro_iso")
cmd.set("stick_transparency", 0.0, "resn CRO")
cmd.disable("tyr143_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 143")
cmd.disable("his146_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 146")
cmd.disable("thr201_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 201")
cmd.disable("ser203_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 203")
cmd.disable("glu220_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 220")
cmd.set_view(
    """
(\
    -0.841515303,   -0.538535833,   -0.042716075,\
    -0.519398630,    0.828276694,   -0.210098386,\
     0.148536488,   -0.154625833,   -0.976732075,\
     0.001062250,   -0.001245134,  -28.732223511,\
    32.638053894,   34.972019196,   38.733093262,\
  -149.667205811,  207.347778320,  -20.000000000 )
"""
)
cmd.png("cro-n.png", dpi=1000)

cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.enable("tyr143_iso")
cmd.set("stick_transparency", 0.0, "resi 143")
cmd.disable("his146_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 146")
cmd.disable("thr201_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 201")
cmd.disable("ser203_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 203")
cmd.disable("glu220_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 220")
cmd.set_view(
    """
(\
    -0.699214101,   -0.229481816,   -0.677073538,\
     0.107669927,    0.902456880,   -0.417058796,\
     0.706748486,   -0.364527941,   -0.606302917,\
    -0.000772315,   -0.001222866,  -22.819011688,\
    38.325408936,   39.157699585,   33.325714111,\
  -155.361190796,  201.653793335,  -20.000000000 )
"""
)
cmd.png("tyr145-n.png", dpi=1000)

cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.disable("tyr143_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 143")
cmd.enable("his146_iso")
cmd.set("stick_transparency", 0.0, "resi 146")
cmd.disable("thr201_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 201")
cmd.disable("ser203_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 203")
cmd.disable("glu220_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 220")
cmd.set_view(
    """
(\
     0.574742734,   -0.673360229,    0.465024173,\
     0.108903490,    0.626136899,    0.772043109,\
    -0.811044395,   -0.393093526,    0.433194339,\
     0.000702321,   -0.000347649,  -27.929138184,\
    39.389068604,   38.298980713,   40.739974976,\
  -150.604400635,  206.410583496,  -20.000000000 )
"""
)
cmd.png("his148-n.png", dpi=1000)


cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.disable("tyr143_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 143")
cmd.disable("his146_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 146")
cmd.enable("thr201_iso")
cmd.set("stick_transparency", 0.0, "resi 201")
cmd.disable("ser203_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 203")
cmd.disable("glu220_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 220")
cmd.set_view(
    """
(\
    -0.358965665,   -0.570977032,    0.738322377,\
    -0.634493172,    0.729414344,    0.255595744,\
    -0.684483767,   -0.376721680,   -0.624129891,\
     0.002785400,   -0.000746295,  -25.533874512,\
    38.447463989,   35.923187256,   39.743259430,\
  -152.961395264,  204.053588867,  -20.000000000 )
"""
)
cmd.png("thr203-n.png", dpi=1000)

cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.disable("tyr143_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 143")
cmd.disable("his146_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 146")
cmd.disable("thr201_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 201")
cmd.enable("ser203_iso")
cmd.set("stick_transparency", 0.0, "resi 203")
cmd.disable("glu220_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 220")
cmd.set_view(
    """
(\
    -0.469252765,   -0.827867150,   -0.307288110,\
    -0.347087085,    0.492887676,   -0.797838151,\
     0.811979473,   -0.267739385,   -0.518630445,\
     0.000395618,    0.000453047,  -21.777130127,\
    39.106773376,   34.605056763,   34.570110321,\
  -156.521820068,  200.493164062,  -20.000000000 )
"""
)
cmd.png("ser205-n.png", dpi=1000)

cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.disable("tyr143_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 143")
cmd.disable("his146_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 146")
cmd.disable("thr201_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 201")
cmd.disable("ser203_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resi 203")
cmd.enable("glu220_iso")
cmd.set("stick_transparency", 0.0, "resi 220")
cmd.set_view(
    """
(\
     0.394451320,   -0.690276980,   -0.606560528,\
     0.803602815,    0.579240561,   -0.136586130,\
     0.445624709,   -0.433566839,    0.783204556,\
    -0.002626797,   -0.000129406,  -33.662319183,\
    41.208389282,   33.075607300,   34.572891235,\
  -144.639602661,  212.375381470,  -20.000000000 )
"""
)
cmd.png("glu222-n.png", dpi=1000)


cmd.enable("cro_iso")
cmd.set("stick_transparency", 0.0, "resn CRO")
cmd.enable("tyr143_iso")
cmd.set("stick_transparency", 0.0, "resi 143")
cmd.enable("his146_iso")
cmd.set("stick_transparency", 0.0, "resi 146")
cmd.enable("thr201_iso")
cmd.set("stick_transparency", 0.0, "resi 201")
cmd.enable("ser203_iso")
cmd.set("stick_transparency", 0.0, "resi 203")
cmd.enable("glu220_iso")
cmd.set("stick_transparency", 0.0, "resi 220")
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
