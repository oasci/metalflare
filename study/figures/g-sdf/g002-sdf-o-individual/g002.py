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
        "cro": "../../../analysis/001-rogfp-md/data/sdf/resnamecro_o.dx",
        "wat": "../../../analysis/001-rogfp-md/data/sdf/watpwire_o.dx",
        "tyr143": "../../../analysis/001-rogfp-md/data/sdf/resid143_o.dx",
        "his146": "../../../analysis/001-rogfp-md/data/sdf/resid146_o.dx",
        "thr201": "../../../analysis/001-rogfp-md/data/sdf/resid201_o.dx",
        "ser203": "../../../analysis/001-rogfp-md/data/sdf/resid203_o.dx",
        "glu220": "../../../analysis/001-rogfp-md/data/sdf/resid220_o.dx",
    },
    "oxd": {
        "cro": "../../../analysis/004-rogfp-oxd-md/data/sdf/resnamecro_o.dx",
        "wat": "../../../analysis/004-rogfp-oxd-md/data/sdf/watpwire_o.dx",
        "tyr143": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid143_o.dx",
        "his146": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid146_o.dx",
        "thr201": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid201_o.dx",
        "ser203": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid203_o.dx",
        "glu220": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid220_o.dx",
    },
    "cu": {
        "cro": "../../../analysis/003-rogfp-cu-md/data/sdf/resnamecro_o.dx",
        "wat": "../../../analysis/003-rogfp-cu-md/data/sdf/watpwire_o.dx",
        "tyr143": "../../../analysis/003-rogfp-cu-md/data/sdf/resid143_o.dx",
        "his146": "../../../analysis/003-rogfp-cu-md/data/sdf/resid146_o.dx",
        "thr201": "../../../analysis/003-rogfp-cu-md/data/sdf/resid201_o.dx",
        "ser203": "../../../analysis/003-rogfp-cu-md/data/sdf/resid203_o.dx",
        "glu220": "../../../analysis/003-rogfp-cu-md/data/sdf/resid220_o.dx",
    },
}

iso_levels = {
    "cro": 0.1,
    "wat": 0.04,
    "tyr143": 0.1,
    "his146": 0.1,
    "thr201": 0.1,
    "ser203": 0.05,
    "glu220": 0.1,
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


TRANSPARENCY_VALUE = 0.5

cmd.enable("cro_iso")
cmd.set("stick_transparency", 0.0, "resn CRO")
cmd.disable("wat_iso")
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
     0.680121481,   -0.356758982,   -0.640427113,\
     0.701607525,    0.570002258,    0.427567095,\
     0.212494507,   -0.740139902,    0.637975276,\
     0.000000000,    0.000000000,  -13.174731255,\
    37.203998566,   38.641998291,   38.016998291,\
  -165.332687378,  191.682296753,  -20.000000000 )
"""
)
cmd.png("cro-oh.png", dpi=1000)
cmd.set_view(
    """
(\
    -0.703570902,   -0.577237546,    0.414460093,\
    -0.652647495,    0.755601943,   -0.055547722,\
    -0.281092286,   -0.309590697,   -0.908359528,\
     0.001810438,   -0.000195116,  -16.810333252,\
    32.588775635,   32.752437592,   36.772533417,\
  -161.677459717,  195.337524414,  -20.000000000 )
"""
)
cmd.png("cro-og1.png", dpi=1000)

cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.enable("wat_iso")
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
    -0.669074118,   -0.100586735,   -0.736352146,\
    -0.044248376,    0.994415581,   -0.095635235,\
     0.741867363,   -0.031418551,   -0.669791639,\
    -0.000842609,   -0.001617131,  -19.710124969,\
    38.974948883,   37.120910645,   37.409934998,\
  -158.781005859,  198.233978271,  -20.000000000 )
"""
)
cmd.png("water-o.png", dpi=1000)

cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.disable("wat_iso")
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
     0.730647206,   -0.682261944,    0.025039384,\
     0.213718995,    0.263394505,    0.940687954,\
    -0.648414493,   -0.681981027,    0.338261306,\
    -0.000361615,    0.001820015,  -16.445045471,\
    36.916408539,   39.897277832,   36.684154510,\
  -162.254104614,  194.760879517,  -20.000000000 )
"""
)
cmd.png("tyr145-o.png", dpi=1000)

cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.disable("wat_iso")
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
     0.274409592,   -0.395127267,   -0.876659632,\
     0.782247603,    0.621917963,   -0.035458703,\
     0.559227824,   -0.676057220,    0.479764134,\
    -0.003097588,    0.001039370,  -17.706823349,\
    40.317687988,   37.705181122,   41.446556091,\
  -160.800720215,  196.214263916,  -20.000000000 )
"""
)
cmd.png("his148-o.png", dpi=1000)


cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.disable("wat_iso")
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
     0.230216607,   -0.464410692,   -0.855169117,\
     0.779711068,    0.613820255,   -0.123443455,\
     0.582243562,   -0.638378382,    0.503428936,\
    -0.002373266,    0.000178115,  -14.972935677,\
    38.473182678,   35.706871033,   38.821453094,\
  -163.628616333,  193.386367798,  -20.000000000 )
"""
)
cmd.png("thr203-o.png", dpi=1000)

cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.disable("wat_iso")
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
    -0.486048996,   -0.638229430,   -0.596976221,\
    -0.316554338,    0.765282452,   -0.560441434,\
     0.814570308,   -0.083434001,   -0.574006259,\
    -0.000057908,   -0.000161140,  -14.880003929,\
    38.682807922,   34.569961548,   35.343654633,\
  -163.049911499,  193.965072632,  -20.000000000 )
"""
)
cmd.png("ser205-o.png", dpi=1000)

cmd.disable("cro_iso")
cmd.set("stick_transparency", TRANSPARENCY_VALUE, "resn CRO")
cmd.disable("wat_iso")
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
    -0.015909592,   -0.977704465,    0.209322944,\
    -0.234800056,    0.207142323,    0.949697137,\
    -0.971898139,   -0.034041606,   -0.232878774,\
     0.001718637,    0.000591275,  -12.979618073,\
    39.552799225,   34.273124695,   37.279449463,\
  -165.913833618,  191.101150513,  -20.000000000 )
"""
)
cmd.png("glu222-o.png", dpi=1000)


cmd.space("rgb")
cmd.enable("cro_iso")
cmd.set("stick_transparency", 0.0, "resn CRO")
cmd.enable("wat_iso")
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
# cmd.quit()
