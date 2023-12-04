import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from loguru import logger


def get_box_lengths(positions: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    r"""Compute lengths of box edges"""
    logger.trace("Computing box lengths")
    box_lengths: npt.NDArray[np.float64] = np.max(positions, axis=0) - np.min(
        positions, axis=0
    )
    logger.trace("Box lengths: {}", box_lengths)
    return box_lengths


def get_box_vectors(positions: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    r"""Vectors of box edges.

    Args:
        positions: Atomic cartesian coordinates.

    Returns:
        Box vectors.
    """
    box_lengths = get_box_lengths(positions)

    logger.trace("Computing box vectors")
    box_vectors = np.zeros((3, 3), dtype=np.float64)
    np.fill_diagonal(box_vectors, box_lengths)
    logger.trace("Box vector:\n{}", box_vectors)

    return box_vectors


def get_com(universe: mda.Universe) -> npt.NDArray[np.float64]:
    r"""Compute the center of mass"""
    logger.info("Computing center of mass (com)")
    com: npt.NDArray[np.float64] = universe.atoms.center_of_mass()
    logger.debug("com: {}", com)
    return com


def get_box_volume(
    positions: npt.NDArray[np.float64],
) -> float:
    r"""Volume of box.

    Args:
        positions: Atomic cartesian coordinates.

    Returns:
        Volume of structures in universe.
    """
    box_vectors = get_box_vectors(positions)
    logger.trace("Computing box volume")
    volume: float = np.linalg.det(box_vectors)
    logger.trace("Box volume: {}", volume)
    return volume
