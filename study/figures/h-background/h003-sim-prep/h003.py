import pymol
from pymol import cmd, util

# https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha

IMG_SIZE = (3.5, 4)  # in
DPI = 300  # 72 for quick; 300 for publication


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
pymol.finish_launching()

# Setup colors
cmd.bg_color("white")

# Load files
cmd.load("../../../data/001-rogfp-md/structures/protein/1JC0-final.pdb", "rogfp2-red")
cmd.load(
    "../../../data/004-rogfp-oxd-md/structures/protein/1JC1-final.pdb", "rogfp2-oxd"
)
cmd.load("../../../data/003-rogfp-cu-md/structures/protein/1JC0-Cu.pdb", "rogfp2-cu")

# Make selections
cmd.select("cro", "resname cro")
cmd.select(
    "cys-sensor",
    "resi 145 or resi 202",
)
cmd.deselect()

# Prep protein
cmd.align("rogfp2-oxd", "rogfp2-red")
cmd.align("rogfp2-cu", "rogfp2-red")

# Set view
cmd.center("cro")

# Display settings
util.performance(0)
cmd.space("cmyk")
cmd.set("ray_trace_fog", 0)
cmd.set("depth_cue", 0)
cmd.set("antialias", 4)
cmd.set("hash_max", 300)
cmd.set("ray_trace_mode", 1)
cmd.set("opaque_background", 0)
cmd.set("cartoon_discrete_colors", 0)
cmd.set("ambient", 0.2)

# Styling
cmd.show("sticks", "cys")
cmd.show("sticks", "resname HOH")
cmd.show("spheres", "element Cu")

cmd.set("cartoon_transparency", 0.3)

cmd.set_color("cartoon-color", hex_to_rgb("#FCFCFC"))
cmd.color("cartoon-color", "rep cartoon")
# Helix cartoon
cmd.set_color("helix-color", hex_to_rgb("#FCFCFC"))
cmd.color("helix-color", "ss H and rep cartoon")
# Sheets cartoon
cmd.set_color("sheet-color", hex_to_rgb("#00b4d8"))
cmd.color("sheet-color", "ss S and rep cartoon")

cmd.set_color("carbon-color", hex_to_rgb("#B8B8B8"))
cmd.color("carbon-color", "element C and not rep cartoon")
util.cnc("cys-sensor")
cmd.color("carbon-color", "element C and cys-sensor and rep sticks")

cmd.set_view(
    """
(\
    -0.459687650,   -0.526690960,    0.715038598,\
    -0.737144887,    0.675323486,    0.023539083,\
    -0.495280445,   -0.516265035,   -0.698686719,\
    -0.000001140,   -0.000002161, -107.157188416,\
     1.669989824,   -0.949798286,   -1.838884473,\
    70.435340881,  143.879058838,  -20.000000000 )
"""
)

cmd.rebuild()
cmd.refresh()

cmd.enable("rogfp2-red")
cmd.disable("rogfp2-oxd")
cmd.disable("rogfp2-cu")
cmd.ray(IMG_SIZE[0] * DPI, IMG_SIZE[1] * DPI)
cmd.png("rogfp2-red-sim.png", dpi=DPI)

cmd.disable("rogfp2-red")
cmd.disable("rogfp2-oxd")
cmd.enable("rogfp2-cu")
cmd.ray(IMG_SIZE[0] * DPI, IMG_SIZE[1] * DPI)
cmd.png("rogfp2-cu-sim.png", dpi=DPI)

cmd.disable("rogfp2-red")
cmd.enable("rogfp2-oxd")
cmd.disable("rogfp2-cu")
cmd.ray(IMG_SIZE[0] * DPI, IMG_SIZE[1] * DPI)
cmd.png("rogfp2-oxd-sim.png", dpi=DPI)

cmd.enable("rogfp2-red")
cmd.enable("rogfp2-oxd")
cmd.enable("rogfp2-cu")

cmd.refresh()
cmd.quit()
