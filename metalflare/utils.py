import importlib
import os
import re

import numpy as np
import numpy.typing as npt
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from scipy.signal import argrelextrema, savgol_filter


def get_obj_from_string(import_string: str) -> object:
    """Retrieves an object based on an import string and object name.

    Args:
        import_string: The import string, starting from the root module, containing
        the desired object. This function would be
        `"metalflare.utils.get_obj_from_string"`.

    Returns:
        The object identified by the import string.
    """
    module_name, obj_name = import_string.rsplit(".", 1)
    module = importlib.import_module(module_name)
    obj = getattr(module, obj_name)
    return obj


def exists_in_array(
    a_slice: npt.NDArray[np.float64 | np.float32],
    array: npt.NDArray[np.float64 | np.float32],
    rtol: float = 1e-05,
    atol: float = 1e-08,
) -> bool:
    r"""Check if ``a_slice`` exists in an ``array``.

    Args:
        a_slice: An example slice of ``array``'s first dimension to check.
        array: Array to check.
        rtol: Relative tolerance for `np.isclose`.
        atol: Absolute tolerance for `np.isclose`.

    Returns:
        If `a_slice` is present in `array` within some tolerance.
    """
    ndim = int(array.ndim)
    exists = np.isclose(a_slice, array, rtol=rtol, atol=atol)
    # For each additional dimension after 1 we check the last dimension
    # if they are all true.
    for _ in range(1, ndim):
        exists = exists.all(axis=-1)
    # At the very end, we will have a 1D array. If any are True, then a_slice
    # exists in array
    return bool(exists.any())


def get_extrema(x, y, extrema_order=3, polyorder=2, window_length=4):
    local_maxima_idxs = argrelextrema(
        savgol_filter(y, window_length=window_length, polyorder=polyorder),
        np.greater,
        order=extrema_order,
    )[0]
    local_maxima = y[local_maxima_idxs]
    local_maxima_x = x[local_maxima_idxs]
    local_minima_idxs = argrelextrema(
        savgol_filter(y, window_length=window_length, polyorder=polyorder),
        np.less,
        order=extrema_order,
    )[0]
    local_minima = y[local_minima_idxs]
    local_minima_x = x[local_minima_idxs]
    local_extrema = np.concatenate([local_maxima, local_minima])
    local_extrema_x = np.concatenate([local_maxima_x, local_minima_x])

    sort_idx = np.argsort(local_extrema_x)
    local_extrema = local_extrema[sort_idx]
    local_extrema_x = local_extrema_x[sort_idx]

    return local_extrema_x, local_extrema


def format_feature_name(feature):
    # Split the feature name into parts
    parts = feature.split("-")

    # Function to convert amino acid names
    def convert_aa_name(name):
        aa_dict = {
            "ala": "Ala",
            "cys": "Cys",
            "asp": "Asp",
            "glu": "Glu",
            "phe": "Phe",
            "gly": "Gly",
            "his": "His",
            "ile": "Ile",
            "lys": "Lys",
            "leu": "Leu",
            "met": "Met",
            "asn": "Asn",
            "pro": "Pro",
            "gln": "Gln",
            "arg": "Arg",
            "ser": "Ser",
            "thr": "Thr",
            "val": "Val",
            "trp": "Trp",
            "tyr": "Tyr",
        }
        match = re.match(r"([a-z]+)(\d+)", name)
        if match:
            aa, num = match.groups()
            if aa.lower() != "cro":
                return f"{aa_dict.get(aa, aa.capitalize())}{int(num)+2}"
            else:
                return f"{aa_dict.get(aa, aa.capitalize())}66"
        return name.capitalize()

    # Function to format atom names
    def format_atom(atom):
        return atom.upper()

    # Process each part
    formatted_parts = []
    for part in parts[
        :-1
    ]:  # All parts except the last one (which is 'dihedral' or 'dist')
        atoms = part.split("_")
        formatted_atoms = [convert_aa_name(atoms[0])] + [
            format_atom(a) for a in atoms[1:]
        ]
        formatted_parts.append("-".join(formatted_atoms))

    # Join the parts
    if parts[-1] == "dihedral":
        return f"{' '.join(formatted_parts)} dihedral"
    elif parts[-1] == "dist":
        return f"{' - '.join(formatted_parts)} distance"
    else:
        return " ".join(formatted_parts)


def load_features(file_paths, normalize_distances=False, transform_dihedrals=False):
    df_list = []
    for path in file_paths:
        data = np.load(path)
        base_name = os.path.basename(path).split(".")[0]

        if "dihedral" in base_name:
            if transform_dihedrals:
                data = (1 - np.cos(data)) / 2
            df = pd.DataFrame({base_name: data})
        elif "dist" in base_name:
            if normalize_distances:
                scaler = MinMaxScaler()
                data = scaler.fit_transform(X=data.reshape(-1, 1)).flatten()
            df = pd.DataFrame({base_name: data})

        df_list.append(df)

    combined_df = pd.concat(df_list, axis=1)
    return combined_df
