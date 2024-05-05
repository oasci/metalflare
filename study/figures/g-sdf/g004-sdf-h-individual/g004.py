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
        "cro": "../../../analysis/001-rogfp-md/data/sdf/resnamecro_h.dx",
        "wat": "../../../analysis/001-rogfp-md/data/sdf/watpwire_h.dx",
        "tyr143": "../../../analysis/001-rogfp-md/data/sdf/resid143_h.dx",
        "his146": "../../../analysis/001-rogfp-md/data/sdf/resid146_h.dx",
        "thr201": "../../../analysis/001-rogfp-md/data/sdf/resid201_h.dx",
        "ser203": "../../../analysis/001-rogfp-md/data/sdf/resid203_h.dx",
        "glu220": "../../../analysis/001-rogfp-md/data/sdf/resid220_h.dx",
    },
    "oxd": {
        "cro": "../../../analysis/004-rogfp-oxd-md/data/sdf/resnamecro_h.dx",
        "wat": "../../../analysis/004-rogfp-oxd-md/data/sdf/watpwire_h.dx",
        "tyr143": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid143_h.dx",
        "his146": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid146_h.dx",
        "thr201": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid201_h.dx",
        "ser203": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid203_h.dx",
        "glu220": "../../../analysis/004-rogfp-oxd-md/data/sdf/resid220_h.dx",
    },
    "cu": {
        "cro": "../../../analysis/003-rogfp-cu-md/data/sdf/resnamecro_h.dx",
        "wat": "../../../analysis/003-rogfp-cu-md/data/sdf/watpwire_h.dx",
        "tyr143": "../../../analysis/003-rogfp-cu-md/data/sdf/resid143_h.dx",
        "his146": "../../../analysis/003-rogfp-cu-md/data/sdf/resid146_h.dx",
        "thr201": "../../../analysis/003-rogfp-cu-md/data/sdf/resid201_h.dx",
        "ser203": "../../../analysis/003-rogfp-cu-md/data/sdf/resid203_h.dx",
        "glu220": "../../../analysis/003-rogfp-cu-md/data/sdf/resid220_h.dx",
    },
}

iso_levels = {
    "cro": 0.6,
    "wat": 0.08,
    "tyr143": 0.25,
    "his146": 0.22,
    "thr201": 0.3,
    "ser203": 0.25,
    "glu220": 0.25,
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
    -0.807500005,   -0.540036440,    0.237262189,\
    -0.531474948,    0.840588927,    0.104459561,\
    -0.255844355,   -0.041757762,   -0.965800405,\
     0.001644215,   -0.001051879,  -26.377836227,\
    32.796859741,   34.064193726,   41.006332397,\
  -151.677536011,  205.337448120,  -20.000000000 )
"""
)
cmd.png("cro-h.png", dpi=1000)

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
    -0.748222709,    0.065092131,   -0.660238981,\
     0.195873305,    0.972470760,   -0.126103342,\
     0.633858442,   -0.223690987,   -0.740377247,\
    -0.001481478,   -0.001338036,  -16.759716034,\
    39.032268524,   37.349628448,   37.025672913,\
  -161.658065796,  195.356918335,  -20.000000000 )
"""
)
cmd.png("water-h.png", dpi=1000)

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
    -0.546161175,   -0.327469409,   -0.771012664,\
     0.286017507,    0.792188942,   -0.539070845,\
     0.787321329,   -0.514957190,   -0.338993728,\
    -0.001563277,   -0.000403425,  -26.318460464,\
    38.368324280,   39.629146576,   34.004287720,\
  -152.129653931,  204.885330200,  -20.000000000 )
"""
)
cmd.png("tyr145-h.png", dpi=1000)

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
     0.363001913,    0.076966047,   -0.928596377,\
     0.843315959,    0.396674275,    0.362540334,\
     0.396244466,   -0.914716125,    0.079084918,\
    -0.002491723,    0.000915487,  -28.869071960,\
    40.539337158,   40.874412537,   38.440757751,\
  -149.647964478,  207.367019653,  -20.000000000 )
"""
)
cmd.png("his148-h.png", dpi=1000)


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
     0.093285978,   -0.613842309,   -0.783887506,\
     0.553596973,    0.686364293,   -0.471596241,\
     0.827518642,   -0.389974266,    0.403861821,\
    -0.002239421,    0.000102347,  -24.747089386,\
    38.710029602,   35.662704468,   39.956199646,\
  -153.748245239,  203.266738892,  -20.000000000 )
"""
)
cmd.png("thr203-h.png", dpi=1000)

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
    -0.396196812,   -0.205815926,   -0.894786000,\
    -0.270139217,    0.957521558,   -0.100631043,\
     0.877504468,    0.201836482,   -0.434967756,\
    -0.000851233,   -0.002343269,  -10.561386108,\
    38.676235199,   34.276218414,   36.444053650,\
  -166.611755371,  190.403228760,  -20.000000000 )
"""
)
cmd.png("ser205-h.png", dpi=1000)

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
    -0.391686231,   -0.380866826,   -0.837560177,\
     0.489106327,    0.684823871,   -0.540145397,\
     0.779306054,   -0.621237397,   -0.081941120,\
    -0.001867004,    0.000305947,  -26.014509201,\
    39.617992401,   31.804334641,   36.105281830,\
  -152.429306030,  204.585678101,  -20.000000000 )
"""
)
cmd.png("glu222-h.png", dpi=1000)


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
     0.601788819,   -0.276795864,   -0.749145329,\
     0.728802443,    0.573926032,    0.373390585,\
     0.326590866,   -0.770695567,    0.547114789,\
    -0.002086747,    0.000230547,  -47.225772858,\
    36.249961853,   36.609130859,   38.500560760,\
  -131.188705444,  225.826278687,  -20.000000000 )
"""
)
cmd.quit()
# cmd.refresh()
