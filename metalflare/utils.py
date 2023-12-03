import importlib

import numpy as np
import numpy.typing as npt


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
