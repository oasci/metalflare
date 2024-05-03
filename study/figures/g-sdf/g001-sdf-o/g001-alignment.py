import os
from collections.abc import Mapping

import numpy as np
import numpy.typing as npt
import pymol
from pymol import cmd
from scipy.optimize import minimize

# https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha

# Start PyMOL session
pymol.finish_launching()


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
    cmd.translate(translate.tolist(), object=mobile_object, camera=0)
    cmd.rotate("x", angles[0], object=mobile_object, camera=0)
    cmd.rotate("y", angles[1], object=mobile_object, camera=0)
    cmd.rotate("z", angles[2], object=mobile_object, camera=0)


def calculate_rmsd(
    x: Mapping[float],
    mobile_object: str,
    ref_object: str,
    atom_selection: str = "resn CRO",
    temp_object_name: str = "temp-moved",
) -> float:
    """Computes RMSD between two PyMOL selections resulting from a proposed spatial
    transformation. This functions copies the `mobile_object` and performs a
    transformation specified by `x` on the temporary object. This keeps the position
    of `mobile_object` unchanged.

    This can be used to optimize the spatial transformation RMSD.

    Args:
        x: Custom transformation vector of
            `[tran_x, trans_y, trans_z, rot_x, rot_y, rot_z]`. This is used in
            `move_object`.
        mobile_object: PyMOL object name to be transformed.
        ref_object: PyMOL object name of the reference object.
        atom_selection: PyMOL atom selection string to compute the RMSD between
            `mobile_object` and `ref_object`.
        temp_object_name: PyMOL object name to call the temporary copy. This will be
            deleted after making the trial move.

    Returns:
        Root-mean-squared deviation of `atom_selection` between `mobile_object` and
        `ref_object`.
    """
    cmd.copy(temp_object_name, mobile_object)
    move_object(temp_object_name, x)
    rmsd = cmd.rms_cur(
        f"{temp_object_name} and ({atom_selection})",
        f"{ref_object} and ({atom_selection})",
    )
    cmd.delete(temp_object_name)
    return rmsd


def opt_object_alignment(
    mobile_object: str, ref_object: str
) -> npt.NDArray[np.floating]:
    """Computes the optimal translation and rotation vector to minimize the RMSD between
    two PyMOL objects.

    Args:
        mobile_object: PyMOL object name to be transformed.
        ref_object: PyMOL object name of the reference object.

    Returns:
        Custom transformation vector of
        `[tran_x, trans_y, trans_z, rot_x, rot_y, rot_z]`.
    """
    # Performs minimization of object translation and rotation
    initial_guess = [0, 0, 0, 0, 0, 0]
    result = minimize(
        calculate_rmsd,
        initial_guess,
        args=(mobile_object, ref_object),
        method="Nelder-Mead",
        bounds=(
            (-30, 30),
            (-30, 30),
            (-30, 30),
            (-180, 180),
            (-180, 180),
            (-180, 180),
        ),
        tol=1e-12,
    )
    move_object(mobile_object, result.x)
    return result.x


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

# We perform two optimizations of each object to the reduced protein.
# The first optimization gets an okay alignment, but we found that performing another
# optimization can fine tune the alignment.
oxd_transform = opt_object_alignment(
    mobile_object="oxd_protein", ref_object="red_protein"
)
oxd_transform = opt_object_alignment(
    mobile_object="oxd_protein", ref_object="red_protein"
)
np.save("transform_oxd.npy", oxd_transform)

cu_transform = opt_object_alignment(
    mobile_object="cu_protein", ref_object="red_protein"
)
cu_transform = opt_object_alignment(
    mobile_object="cu_protein", ref_object="red_protein"
)
np.save("transform_cu.npy", cu_transform)

cmd.quit()
