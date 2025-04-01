import os

import pymol
from pymol import cmd, util

# https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha

os.chdir(os.path.dirname(os.path.realpath(__file__)))

DIR_BASE = "../../../"


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

sele_relevant = "(resn CRO or resi 201 or resi 203 or resi 220)"
for i in range(0, 11):
    # Load structure
    label_cluster = chr(i + 65)
    path_pdb = os.path.join(
        DIR_BASE,
        f"analysis/009-pw-configs/data/cluster_pdbs/cluster_{label_cluster}.pdb"
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
    cmd.remove(
        f"byres resn WAT and not (resn WAT within 8 of resn CRO) and model cluster_{label_cluster}"
    )

    sele_mobile = sele_relevant + f" and model cluster_{label_cluster}"
    cmd.show(representation="licorice", selection=sele_mobile)

    # Align relevant residues
    if i != 0:
        sele_target = sele_relevant + " and model cluster_A"
        cmd.align(sele_mobile, sele_target)


# Hide nonpolar residues
cmd.hide("h. and (e. c extend 1)")
cmd.center("resi 203 and model cluster_A")

# pymol.refresh()

