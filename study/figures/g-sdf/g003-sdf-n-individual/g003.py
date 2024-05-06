from collections.abc import Mapping

import numpy as np
import pymol
from pymol import cmd, util

# https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha


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

iso_levels = {
    "cro": 0.09,
    "tyr143": 0.09,
    "his146": 0.09,
    "thr201": 0.09,
    "ser203": 0.09,
    "glu220": 0.09,
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
        cmd.isomesh(name=iso_label, map=map_label)
        cmd.isolevel(iso_label, iso_levels[res_label])
        cmd.group(name=f"{res_label}_iso", members=iso_label, action="add")

        cmd.color(f"{sys_label}_color", iso_label)
        cmd.set("transparency", 0.2, iso_label)


# Disable other objects
cmd.disable("oxd_protein")
cmd.disable("cu_protein")
cmd.disable("maps")

# Set view
cmd.color("grey70", "element C")

# Figure settings
util.performance(0)
cmd.space("cmyk")
cmd.set("mesh_width", 1.5)
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
    -0.713077724,   -0.701010048,    0.010032250,\
    -0.626331806,    0.643409133,    0.440112710,\
    -0.314971328,    0.307547778,   -0.897879541,\
     0.001500354,   -0.001301490,  -15.251414299,\
    30.964843750,   33.253074646,   38.581691742,\
  -163.248641968,  193.766342163,  -20.000000000 )
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
    -0.775559962,   -0.575533450,    0.259353071,\
    -0.586101055,    0.809085667,    0.042806949,\
    -0.234468073,   -0.118816726,   -0.964825809,\
     0.000000000,    0.000000000,  -15.039458275,\
    41.667057037,   39.947574615,   29.875352859,\
  -163.468048096,  193.546936035,  -20.000000000 )
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
     0.338627607,    0.308200926,   -0.889008164,\
     0.090006888,    0.929875195,    0.356649607,\
     0.936592400,   -0.200803891,    0.287140876,\
    -0.001674286,   -0.001392300,  -22.638912201,\
    39.828647614,   41.086414337,   36.963531494,\
  -155.992446899,  201.022537231,  -20.000000000 )
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
     0.858347416,    0.171213508,   -0.483652413,\
    -0.046125349,    0.964596391,    0.259603322,\
     0.510979593,   -0.200536281,    0.835861146,\
    -0.000913270,   -0.001376533,  -17.563423157,\
    40.380111694,   35.601814270,   41.503143311,\
  -160.975646973,  196.039337158,  -20.000000000 )
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
     0.755386353,   -0.059837028,    0.652534127,\
    -0.642282307,    0.129595384,    0.755410016,\
    -0.129761919,   -0.989751458,    0.059450448,\
     0.000751109,    0.001053762,  -16.790933609,\
    40.718509674,   34.968013763,   35.085838318,\
  -161.790435791,  195.224578857,  -20.000000000 )
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
     0.086835407,   -0.498051465,    0.862782896,\
    -0.941057384,    0.243140176,    0.235075340,\
    -0.326848060,   -0.832350969,   -0.447593719,\
     0.000049660,    0.000069670,  -18.572540283,\
    39.984493256,   30.044538498,   33.884197235,\
  -159.939743042,  197.075271606,  -20.000000000 )
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
