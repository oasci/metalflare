import pymol
from pymol import cmd, util

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
    return (r / 255.0, g / 255.0, b / 255.0)


# Start PyMOL session
pymol.finish_launching()

# Setup colors
cmd.bg_color("white")

# Load files
cmd.load("https://files.rcsb.org/download/2Y0G.pdb", "egfp")
cmd.load("https://files.rcsb.org/download/1JC0.pdb", "rogfp2-reduced")
cmd.load("https://files.rcsb.org/download/1JC1.pdb", "rogfp2-oxidized")
cmd.load("../../../data/003-rogfp-cu-md/structures/protein/1JC0-Cu.pdb", "rogfp2-cu")

# Make selections
cmd.remove("not chain A and not chain X")
cmd.select("cro", "resname cro")
cmd.select(
    "cys-sensor",
    "((resi 147 or resi 204) and not model rogfp2-cu) or ((resi 145 or resi 202) and model rogfp2-cu)",
)
cmd.select(
    "cro-water",
    "(model rogfp2-reduced and resi 266) or (model rogfp2-oxidized and resi 241) or (model egfp and resi 2084)",
)
cmd.deselect()

# Remove unnecessary waters
cmd.remove("(resname HOH or resname WAT) and not cro-water")

# Prep protein
cmd.h_add("not model rogfp2-cu")
cmd.remove("resname CRO and name H15")
cmd.remove("element H and model rogfp2-cu and cys-sensor")
cmd.align("rogfp2-oxidized", "rogfp2-reduced")
cmd.align("egfp", "rogfp2-reduced")
cmd.align("rogfp2-cu", "rogfp2-reduced")

# Set view
cmd.center("cys-sensor")
cmd.set_view(
    """
(\
     0.427638769,    0.786577761,    0.445441961,\
     0.903749168,   -0.361738890,   -0.228857800,\
    -0.018881869,    0.500437379,   -0.865563154,\
    -0.000326671,   -0.000341356,  -46.396766663,\
   185.186950684,    8.843163490,   43.544662476,\
    -8.412878990,  101.244674683,  -20.000000000 )
"""
)

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
cmd.set("ambient", 0)

# Styling
cmd.show("sticks", "cys")
cmd.show("sticks", "resname HOH")
cmd.show("spheres", "element Cu")

cmd.set("cartoon_transparency", 0.8)

cmd.set_color("cartoon-color", hex_to_rgb("#F2F2F2"))
cmd.color("cartoon-color", "rep cartoon")
# Helix cartoon
cmd.set_color("helix-color", hex_to_rgb("#F2F2F2"))
cmd.color("helix-color", "ss H and rep cartoon")
# Sheets cartoon
cmd.set_color("sheet-color", hex_to_rgb("#00b4d8"))
cmd.color("sheet-color", "ss S and rep cartoon")

cmd.set_color("carbon-color", hex_to_rgb("#B8B8B8"))
cmd.color("carbon-color", "element C and not rep cartoon")
cmd.color("carbon-color", "element C and cys-sensor")

util.cnc("cys-sensor")

cmd.rebuild()
cmd.refresh()

# Pictures
cmd.enable("egfp")
cmd.disable("rogfp2-reduced")
cmd.disable("rogfp2-oxidized")
cmd.disable("rogfp2-cu")
cmd.png("egfp.png", dpi=1000)
cmd.disable("egfp")
cmd.enable("rogfp2-reduced")
cmd.disable("rogfp2-oxidized")
cmd.disable("rogfp2-cu")
cmd.png("rogfp2-reduced.png", dpi=1000)
cmd.disable("egfp")
cmd.disable("rogfp2-reduced")
cmd.enable("rogfp2-oxidized")
cmd.disable("rogfp2-cu")
cmd.png("rogfp2-oxidized.png", dpi=1000)

cmd.disable("egfp")
cmd.disable("rogfp2-reduced")
cmd.disable("rogfp2-oxidized")
cmd.enable("rogfp2-cu")
cmd.png("rogfp2-cu.png", dpi=1000)

cmd.enable("rogfp2-reduced")
cmd.enable("rogfp2-oxidized")
cmd.disable("rogfp2-cu")

cmd.refresh()
