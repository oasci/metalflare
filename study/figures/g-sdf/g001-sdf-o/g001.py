import pymol
from pymol import cmd, stored
from scipy.optimize import minimize

# https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha

ISO_VALUE = 0.0018

# Start PyMOL session
pymol.finish_launching()


def calculate_rmsd(x, mobile_object, ref_object):
    cmd.copy("trans-move", mobile_object)
    move_object("trans-move", x)
    rmsd = cmd.rms_cur(f"trans-move and resn CRO", f"{ref_object} and resn CRO")
    cmd.delete("trans-move")
    return rmsd


def move_object(mobile_object, x):
    translate = x[:3]
    angles = x[3:]
    cmd.translate(translate.tolist(), object=mobile_object, camera=0)
    cmd.rotate("x", angles[0], object=mobile_object, camera=0)
    cmd.rotate("y", angles[1], object=mobile_object, camera=0)
    cmd.rotate("z", angles[2], object=mobile_object, camera=0)


def calc_object_matrix(mobile_object, ref_object):
    # Performs minimization of object translation and rotation
    initial_guess = [0, 0, 0, 0, 0, 0]
    result = minimize(
        calculate_rmsd,
        initial_guess,
        args=(mobile_object, ref_object),
        method="Nelder-Mead",
        bounds=(
            (-50, 50),
            (-50, 50),
            (-50, 50),
            (-360, 360),
            (-360, 360),
            (-360, 360),
        ),
        tol=1e-12,
    )
    move_object(mobile_object, result.x)
    return result.x


# Setup colors
cmd.bg_color("white")
cmd.color("grey70", "element C")
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

cmd.select("red_cro", "red_protein and resn CRO")
cmd.select("oxd_cro", "oxd_protein and resn CRO")
cmd.select("cu_cro", "cu_protein and resn CRO")

# Display specific residues as sticks
cmd.select("proton_wire", "(resi 143 or resi 146 or resi 201 or resi 220)")
cmd.show("sticks", "proton_wire")
cmd.select("cys", "(resi 145 or resi 202)")
cmd.show("sticks", "cys")
cmd.hide("cartoon", "(not proton_wire and not cys)")

# Select CRO residues
cmd.select("cro_red", "(resn CRO and model red_protein)")
cmd.select("cro_oxd", "(resn CRO and model oxd_protein)")
cmd.select("cro_cu", "(resn CRO and model cu_protein)")

# Load density maps
cmd.load("../../../analysis/001-rogfp-md/data/sdf/resid65_oh-o.dx", "red")
cmd.load("../../../analysis/004-rogfp-oxd-md/data/sdf/resid65_oh-o.dx", "oxd")
cmd.load("../../../analysis/003-rogfp-cu-md/data/sdf/resid65_oh-o.dx", "cu")

# Generate and adjust surfaces
cmd.isomesh("red_iso", "red", 1.0)
cmd.isomesh("oxd_iso", "oxd", 1.0)
cmd.isomesh("cu_iso", "cu", 1.0)

cmd.isolevel("red_iso", ISO_VALUE)
cmd.isolevel("oxd_iso", ISO_VALUE)
cmd.isolevel("cu_iso", ISO_VALUE)

# Change colors of isosurfaces
cmd.color("red_color", "red_iso")
cmd.color("oxd_color", "oxd_iso")
cmd.color("cu_color", "cu_iso")

# Perform the alignment using the align command which returns the RMSD and the transformation matrix
oxd_transform = calc_object_matrix(
    mobile_object="oxd_protein", ref_object="red_protein"
)
cu_transform = calc_object_matrix(mobile_object="cu_protein", ref_object="red_protein")
move_object("oxd_iso", oxd_transform)
move_object("cu_iso", cu_transform)

oxd_transform = calc_object_matrix(
    mobile_object="oxd_protein", ref_object="red_protein"
)
cu_transform = calc_object_matrix(mobile_object="cu_protein", ref_object="red_protein")
move_object("oxd_iso", oxd_transform)
move_object("cu_iso", cu_transform)

# Disable other proteins
# cmd.disable("oxd_protein")
# cmd.disable("cu_protein")

# Ensure PyMOL updates the display
cmd.refresh()
