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
cmd.set_color("cartoon-color", hex_to_rgb("#F2F2F2"))
cmd.set_color("helix-color", hex_to_rgb("#F2F2F2"))
cmd.set_color("sheet-color", hex_to_rgb("#00CCF5"))
cmd.set_color("helix-color", hex_to_rgb("#FFE45E"))
cmd.set_color("loop-color", hex_to_rgb("#F2F2F2"))
cmd.set_color("carbon-color", hex_to_rgb("#B8B8B8"))

cmd.set_color("reduced-color", hex_to_rgb("#1e2e79"))
cmd.set_color("oxidized-color", hex_to_rgb("#EC4067"))
cmd.set_color("cu-color", hex_to_rgb("#f99752"))

# Load files
cmd.load(
    "../../../figures/b-cys/b004-cys147_ca-cys204_ca/pdbs/f001-reduced-4.310.pdb",
    "reduced",
)
cmd.load(
    "../../../figures/b-cys/b004-cys147_ca-cys204_ca/pdbs/f001-oxidized-4.080.pdb",
    "oxidized",
)
cmd.load("../../../figures/b-cys/b004-cys147_ca-cys204_ca/pdbs/f001-cu-4.310.pdb", "cu")

# Make selections
residue_selection = "(resi 143-146 or resi 201-203 or resi 220)"
cmd.deselect()
cmd.select("cro", "resname cro")
cmd.select(
    "relevant-reduced",
    f"{residue_selection} and model reduced",
)
cmd.create(name="relevant-reduced", selection="relevant-reduced")
cmd.select(
    "relevant-oxidized",
    f"{residue_selection} and model oxidized",
)
cmd.create(name="relevant-oxidized", selection="relevant-oxidized")
cmd.select(
    "relevant-cu",
    f"{residue_selection} and model cu",
)
cmd.create(name="relevant-cu", selection="relevant-cu")
cmd.deselect()

# Set view
cmd.center("relevant-reduced")
cmd.set_view(
    """
(-0.3279060423374176, -0.20367984473705292, 0.9224839806556702, -0.3389420807361603, 0.9368278384208679, 0.08636336028575897, -0.8818066716194153, -0.2843515872955322, -0.37623244524002075, 0.0007185516878962517, -0.00032091327011585236, -47.9632453918457, 35.017723083496094, 38.951385498046875, 30.955150604248047, 35.858924865722656, 59.94841003417969, -20.0)
"""
)

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
cmd.set("label_size", 24)
cmd.set("label_font_id", 7)

# Styling
cmd.show("sticks", "relevant-reduced")
cmd.show("sticks", "relevant-oxidized")
cmd.show("sticks", "relevant-cu")
# cmd.show("spheres", "element Cu")
cmd.hide(representation="sticks", selection="element H")
cmd.show(representation="sticks", selection="(resi 145 or resi 202) and name HG")
cmd.show(representation="sticks", selection="(resi 143) and name HH")
cmd.show(representation="sticks", selection="(resi 146) and name HD1")
cmd.show(representation="sticks", selection="(resi 201) and name HG1")
cmd.show(representation="sticks", selection="(resi 203) and name HG")
cmd.show(representation="sticks", selection="(resi 220) and name HE2")
cmd.hide(representation="cartoon", selection="relevant-reduced")
cmd.hide(representation="cartoon", selection="relevant-oxidized")
cmd.hide(representation="cartoon", selection="relevant-cu")

cmd.set("cartoon_transparency", 0.6)

cmd.color("cartoon-color", "rep cartoon")
cmd.color("helix-color", "ss H and rep cartoon")
cmd.color("sheet-color", "ss S and rep cartoon")
cmd.color("helix-color", "ss h and rep cartoon")
cmd.color("loop-color", "ss l and rep cartoon")
cmd.color("carbon-color", "element C and rep sticks")

util.cnc("relevant-reduced")
util.cnc("relevant-oxidized")
util.cnc("relevant-cu")

cmd.label(
    '''(name CA+C1*+C1' and (byres(relevant-reduced)))''','''"%s-%s"%(resn,resi)'''
)
cmd.label(
    '''(name CA+C1*+C1' and (byres(relevant-oxidized)))''','''"%s-%s"%(resn,resi)'''
)
cmd.label(
    '''(name CA+C1*+C1' and (byres(relevant-cu)))''','''"%s-%s"%(resn,resi)'''
)
cmd.set("label_position", (0, -2, 3))

cmd.rebuild()
cmd.refresh()

# Pictures
cmd.disable("all")
cmd.enable("reduced")
cmd.enable("relevant-reduced")
cmd.png("relevant-reduced.png", dpi=1000)

cmd.disable("all")
cmd.enable("oxidized")
cmd.enable("relevant-oxidized")
cmd.png("relevant-oxidized.png", dpi=1000)

cmd.disable("all")
cmd.enable("cu")
cmd.enable("relevant-cu")
cmd.png("relevant-cu.png", dpi=1000)

cmd.disable("all")
cmd.enable("reduced")
cmd.enable("relevant-reduced")

cmd.refresh()
