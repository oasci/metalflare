import os

import pymol
from pymol import cmd, util

os.chdir(os.path.abspath(os.path.dirname(__file__)))

# https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha


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
    return (r, g, b)

# Start PyMOL session
pymol.finish_launching()

# Setup colors
cmd.bg_color("white")

SAVE_PNG = True

# Load files
cmd.load("../../../analysis/005-rogfp-glh-md/data/traj/frame_89526.pdb", "reduced")
cmd.load("../../../analysis/007-rogfp-oxd-glh-md/data/traj/frame_50564.pdb", "oxidized")
cmd.load("../../../analysis/006-rogfp-cu-glh-md/data/traj/frame_100971.pdb", "cu")
cmd.center("reduced")

# Remove unnecessary waters
cmd.remove("(resname HOH or resname WAT)")
cmd.remove("element Na")
cmd.remove("element Cl")

# Align
cmd.align("model oxidized", "model reduced")
cmd.align("model cu", "model reduced")

# Make selections
cmd.deselect()
cmd.select(
    "cys-sensor-reduced",
    "(resi 145 or resi 202) and model reduced",
)
cmd.create(name="cys-reduced", selection="cys-sensor-reduced")
cmd.select(
    "cys-sensor-oxidized",
    "(resi 145 or resi 202) and model oxidized",
)
cmd.create(name="cys-oxidized", selection="cys-sensor-oxidized")
cmd.select(
    "cys-sensor-cu",
    "(resi 145 or resi 202) and model cu",
)
cmd.create(name="cys-cu", selection="cys-sensor-cu")
cmd.deselect()

# Display settings
util.performance(0)
cmd.space("rgb")
cmd.set("ray_trace_fog", 0)
cmd.set("depth_cue", 0)
cmd.set("antialias", 4)
cmd.set("hash_max", 300)
cmd.set("ray_trace_mode", 0)
cmd.set("opaque_background", 0)
cmd.set("cartoon_discrete_colors", 0)
cmd.set("direct", 0.5)
cmd.set("ambient", 0.25)
cmd.set("specular", 0)
cmd.set("light_count", 4)
cmd.set("ray_shadows", 0)

# Styling
cmd.show("sticks", "cys-reduced")
cmd.show("sticks", "cys-oxidized")
cmd.show("sticks", "cys-cu")
cmd.show("sticks", "resname HOH")
cmd.show("spheres", "element Cu")

# cmd.set("cartoon_transparency", 0.8)

cmd.set_color("cartoon-color", hex_to_rgb("#F2F2F2"))
cmd.color("cartoon-color", "rep cartoon")
# Helix cartoon
cmd.set_color("helix-color", hex_to_rgb("#F2F2F2"))
cmd.color("helix-color", "ss H and rep cartoon")
# Sheets cartoon
cmd.set_color("sheet-color", hex_to_rgb("#00CCF5"))
cmd.color("sheet-color", "ss S and rep cartoon")
cmd.set_color("helix-color", hex_to_rgb("#FFE45E"))
cmd.color("helix-color", "ss h and rep cartoon")
cmd.set_color("loop-color", hex_to_rgb("#F2F2F2"))
cmd.color("loop-color", "ss l and rep cartoon")
# Carbon color
cmd.set_color("carbon-color", hex_to_rgb("#B8B8B8"))
cmd.color("carbon-color", "element C and rep sticks")

util.cnc("cys-reduced")
util.cnc("cys-oxidized")
util.cnc("cys-cu")

cmd.set_view(
    """
(-0.4751238226890564, 0.4551173746585846, 0.7530767917633057, 0.680943489074707, 0.7322223782539368, -0.01290218811482191, -0.5572898983955383, 0.5066716074943542, -0.6578046679496765, 0.00015270337462425232, -0.00016066618263721466, -31.541841506958008, 35.07742691040039, 46.132991790771484, 36.379547119140625, -14.140244483947754, 77.19920349121094, -20.0)
"""
)

if SAVE_PNG:
    # Pictures
    png_path = "cys-sensor-reduced.png"
    print(f"Rendering {png_path}")
    cmd.disable("all")
    cmd.enable("reduced")
    cmd.enable("cys-reduced")
    cmd.png(png_path, dpi=1000)

    png_path = "cys-sensor-oxidized.png"
    print(f"Rendering {png_path}")
    cmd.disable("all")
    cmd.enable("oxidized")
    cmd.enable("cys-oxidized")
    cmd.png(png_path, dpi=1000)

    png_path = "cys-sensor-cu.png"
    print(f"Rendering {png_path}")
    cmd.disable("all")
    cmd.enable("cu")
    cmd.enable("cys-cu")
    cmd.png(png_path, dpi=1000)

    cmd.enable("reduced")
    cmd.enable("oxidized")
    cmd.enable("cu")

