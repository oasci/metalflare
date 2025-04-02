import os

import pymol
from pymol import cmd, util

# https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha

os.chdir(os.path.dirname(os.path.realpath(__file__)))

DIR_BASE = "../../../"

render_figures = False


def hex_to_rgb(hex_color):
    """
    Convert a hexadecimal color code to a tuple of RGB values normalized between 0 and 1 for PyMOL.

    Args:
        hex_color (str): The hexadecimal color code starting with '#'.

    Returns:
        tuple: A tuple of three floats representing the RGB values.
    """
    # Remove the '#' character if it's there
    hex_color = hex_color.strip("#")

    # Convert the hex values to integer
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)

    # Normalize the RGB values to 0-1 range (PyMOL RGB format)
    return (r / 255.0, g / 255.0, b / 255.0)


# Start PyMOL session
if not render_figures:
    pymol.finish_launching()

# Setup colors
cmd.bg_color("white")

# Display settings
util.performance(0)
cmd.space("rgb")
cmd.set("ray_trace_fog", 0)
cmd.set("depth_cue", 1)
cmd.set("antialias", 4)
cmd.set("hash_max", 300)
cmd.set("ray_trace_mode", 0)
cmd.set("opaque_background", 0)
cmd.set("cartoon_discrete_colors", 0)
cmd.set("ambient", 0)

# Load files

sele_relevant = "(resn CRO or resi 201 or resi 203 or resi 220)"
for i in range(0, 11):
    # Load structure
    label_cluster = chr(i + 65)
    path_pdb = os.path.join(
        DIR_BASE,
        f"analysis/009-pw-configs/data/cluster_pdbs/cluster_{label_cluster}.pdb",
    )
    label_object = f"cluster_{label_cluster}"
    cmd.load(path_pdb, label_object)

    # Fix PBC and move to [0, 0, 0]
    com_cro = cmd.centerofmass(selection=f"resn CRO and model {label_object}")
    cmd.pbc_wrap(label_object, center=com_cro)
    cmd.translate(
        vector=[-val for val in com_cro], selection=f"model {label_object}", camera=0
    )

    cmd.remove("resn Na\+ or resn Cl\-")
    # cmd.remove(
    #     f"byres resn WAT and not (resn WAT within 8 of resn CRO) and model cluster_{label_cluster}"
    # )
    cmd.remove("resn WAT")

    sele_mobile = sele_relevant + f" and model cluster_{label_cluster}"
    cmd.show(representation="licorice", selection=sele_mobile)

    # Align relevant residues
    if i != 0:
        sele_target = sele_relevant + " and model cluster_A"
        cmd.align(sele_mobile, sele_target)


# Hide nonpolar residues
cmd.hide("h. and (e. c extend 1)")
cmd.center("resi 203 and model cluster_A")

cmd.hide(
    representation="cartoon",
    selection="resi 130-135 or resi 164-176 or resi 56-62 or resi 137-148",
)
cmd.hide(selection="resn CU1")

cmd.set_view(
    """
(\
    -0.072506048,    0.982876003,    0.169377849,\
     0.608668327,   -0.090929076,    0.788197279,\
     0.790103614,    0.160242677,   -0.591656387,\
     0.000014191,    0.000051068,  -36.410449982,\
    -5.155659199,   -0.965501785,   -2.455146551,\
    28.928140640,   43.900684357,  -20.000000000 )
"""
)

# Color and render clusters
dpi = 600
width = 4  # in
height = 4  # in

cmd.set_color("DEF-color", hex_to_rgb("#5C72D6"))
cmd.color(
    "DEF-color", "element C & (model cluster_D | model cluster_E | model cluster_F)"
)

cmd.set_color("ABC-color", hex_to_rgb("#F492A9"))
cmd.color(
    "ABC-color", "element C & (model cluster_A | model cluster_B | model cluster_C)"
)

cmd.set_color("JK-color", hex_to_rgb("#4BAA9D"))
cmd.color("JK-color", "element C & (model cluster_J | model cluster_K)")

cmd.set_color("HI-color", hex_to_rgb("#FBC49D"))
cmd.color("HI-color", "element C & (model cluster_H | model cluster_I)")

cmd.set_color("G-color", hex_to_rgb("#90BE6D"))
cmd.color("G-color", "element C & (model cluster_G)")

if render_figures:
    cmd.disable("all")
    cmd.enable("cluster_D, cluster_E, cluster_F")
    cmd.ray(width=width * dpi, height=height * dpi)
    cmd.png("clusters-DEF.png", dpi=dpi)

    cmd.disable("all")
    cmd.enable("cluster_A, cluster_B, cluster_C")
    cmd.ray(width=width * dpi, height=height * dpi)
    cmd.png("clusters-ABC.png", dpi=dpi)

    cmd.disable("all")
    cmd.enable("cluster_J, cluster_K")
    cmd.ray(width=width * dpi, height=height * dpi)
    cmd.png("clusters-JK.png", dpi=dpi)

    cmd.disable("all")
    cmd.enable("cluster_H, cluster_I")
    cmd.ray(width=width * dpi, height=height * dpi)
    cmd.png("clusters-HI.png", dpi=dpi)

    cmd.disable("all")
    cmd.enable("cluster_G")
    cmd.ray(width=width * dpi, height=height * dpi)
    cmd.png("clusters-G.png", dpi=dpi)
